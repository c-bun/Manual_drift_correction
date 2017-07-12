package de.mpicbg.scf.ManualDriftCorrection;

/**
 * Author: Benoit Lombardot, Scientific computing facility @ MPI-CBG  
 * date: 2015-06-03
 */


/**

for the manual_registration:
- load the manual registration script in fiji script editor
- load the image to analyze
- save a set of ROI (point, line or polyline) describing your object movement to the ROI manager
  (in menu "Analyze>tools>ROI manager", use the key "t" to add a region to the ROI manager).
  Only one type of roi can be used for a given dataset
- check with the function "Image>hyperstack>stack to hyperstack" that the time dimension of your data
  matches to the time dimension in FIJI
- Set the time slider to the reference image (all other image of the stack will be registered to this one)
- click run in the script editor. The output will be a registered stack.

 rk  : the ROI in the Roi manager should be in ascending time order (smallest time first)
 rk2 : the plugin currently does 2D translation correction. it can very easily be adapted to correct for 
       rigid transformation or similarity transformation
 rk3 : there should be the same number of landmark at each time step
 rk4 : there should be a landmark at first and last time step. if not the landmark appearing first (last) will be used
       on all anterior (posterior) timestep

 */

import ij.plugin.PlugIn;
import ij.IJ;
import ij.ImagePlus;
import ij.process.ImageProcessor;
import ij.plugin.frame.RoiManager;
import ij.gui.Roi;
import ij.gui.Line;
import ij.gui.PointRoi;
import ij.gui.PolygonRoi;
import ij.CompositeImage;
import ij.ImageStack;



import java.util.ArrayList;
import java.util.Arrays;

import mpicbg.models.Point;
import mpicbg.models.PointMatch;
import mpicbg.ij.InverseTransformMapping;
import mpicbg.ij.Mapping;
import mpicbg.models.InverseCoordinateTransform;
import mpicbg.models.TranslationModel2D;
//import mpicbg.models.RigidModel2D;
//import mpicbg.models.SimilarityModel2D;
//import mpicbg.models.AffineModel2D;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.Affine2D;





/*
 * TODO:
 * - create a proper class: to set land marks, interpolate them, set the model, perform the registration
 * - interface : ref slice, slice range, methods, origin stack (for the landmark czt position proper calculus)
 * - register along z, t or channel
 * (- changing number of landmark  by creating some landmark trajectory)
 * (- handle 3D reg) not sure it worth
 */
 
public class Manual_Registration implements PlugIn {

		boolean interpolate;
		boolean showMatrix;
		
		public void run(String arg) {
			interpolate=true;
			showMatrix=false;
		
			IJ.log("-------------Drift Correction Plugin ------------");	
			ImagePlus imp =  IJ.getImage();
			int[] inputDims = imp.getDimensions();
			int slice_per_frame = inputDims[2]*inputDims[3]; // Nchannel*Nz
			IJ.log("slice_per_frame "+slice_per_frame);
			
			int nFrames = inputDims[4];
			
			RoiManager manager = RoiManager.getInstance();
			if (manager == null){
	    		IJ.log("no manager open");
	    		return;
			}
			Roi[] rois = manager.getRoisAsArray();
			if (rois.length == 1 && rois[0] instanceof PointRoi) {
				PointRoi pointRoi = (PointRoi) rois[0];
				if (pointRoi.getCount(0) > 1) {
					IJ.showMessage("Multi-point selections only supported if one selection per frame");
					return;
				}
			}

			ArrayList<ArrayList<Point>> keyLandmark_list = new ArrayList<ArrayList<Point>>();
			ArrayList<Integer> keyFrame_list = new ArrayList<Integer>();
			int npt = 0;
			IJ.log("n rois : "+rois.length);
			for (int i = 0; i < rois.length; i++) 
			{
				Roi roi = manager.getRoi(i);
				String label = roi.getName();
				if (label==null){ IJ.log("label is null");}
				
				int frame;
				if (imp.isHyperStack())
				{	frame = roi.getTPosition();  }
				else
				{	//int frame = (int) (manager.getSliceNumber(label)/slice_per_frame)+1;
					frame = (int) ((float)(roi.getPosition()-1)/slice_per_frame )+1;
				}
	        	IJ.log("rois: name - frame : "+roi.getName()+" - "+frame);
	        	//IJ.log("roimanager: slice - name : "+ manager.getSliceNumber(label) +" - "+label);
	        	
	        	// Read landmark information associated with key-frames
	        	keyFrame_list.add(frame);
	        	if (roi.getType() == Roi.LINE)
	        	{
	        		npt=2;
	        		Line roi2 = (Line) roi;
	        		Point p;
	        		ArrayList<Point> points = new ArrayList<Point>();
	        		p = new Point( new double[] {(double)roi2.x1d,(double)roi2.y1d} );
	        		points.add(p);
	        		p = new Point( new double[] {(double)roi2.x2d,(double)roi2.y2d} );
	        		points.add(p);
	        		keyLandmark_list.add(points);	        		
	        	}
	            else if(roi.getType() == Roi.POINT |  roi.getType() == Roi.POLYLINE)
	        	{
	        		PolygonRoi roi2 = (PolygonRoi) roi;
	        		int npoints  = roi2.getNCoordinates();
	        		
	        		float[] aux = roi.getFloatPolygon().xpoints;
	        		double[] x_points = new double[aux.length];
	        		for( int j=0; j<aux.length; j++) { x_points[j] = (double)aux[j]; } 
	        		
	        		aux = roi.getFloatPolygon().ypoints;
	        		double[] y_points = new double[aux.length];
	        		for( int j=0; j<aux.length; j++) { y_points[j] = (double)aux[j]; } 

	        		Point p;
	        		ArrayList<Point> points = new ArrayList<Point>();
	        		for (int j=0; j<npoints; j++)
	        		{
		        		p = new Point( new double[] {x_points[j],y_points[j]} );
		        		points.add(p);
	        		}
	        		keyLandmark_list.add(points);

	        		npt = npoints;
	        	}	        
			}
			//IJ.log("key landmark " + keyLandmark_list.toString() );
			//IJ.log("key Frame " + keyFrame_list.toString() );



			// todo: sort keyframe in increasing order
			//		 sort keylandmark according to keyframe
			
			// create a linear interpolated list of point
			ArrayList< ArrayList<Point> > Landmark_list = new ArrayList< ArrayList<Point> >(nFrames);
			int t1_idx = 0; //keyFrame_list.size()-1 ;
			int t2_idx = 1; //t1_idx-1;
			int t1_idx_max = keyFrame_list.size()-1;
			int t1,t2;
			
			int nd=2; // number of dimension of the point
			double[][] p1=new double[npt][];
			double[][] p2=new double[npt][];
			for(int i=0; i<npt; i++){
				//p[i] = new float[nd];
				p1[i] = new double[nd];
				p2[i] = new double[nd];
			}

			t1 = keyFrame_list.get(t1_idx).intValue();
			t2 = keyFrame_list.get(t2_idx).intValue();
			
			for(int t=1; t<=nFrames; t++)
			{	
				if (t2<=t)
				{
					t1=t2;
					t1_idx++;
					t2_idx++;
					if (t1_idx<t1_idx_max)
						t2 = keyFrame_list.get(t2_idx).intValue();
					else{
						t2=t1;
						t1_idx=t1_idx_max; t2_idx=t1_idx_max;
					}
				}
				//IJ.log("interp range ; "+t1+" "+t+" "+t2);
				//IJ.log("t1 idx : " + t1_idx);
				
				ArrayList<Point> point_list = new ArrayList<Point>();
				if (t<=t1){
					for(int pt=0; pt<npt; pt++)
						point_list.add( new Point( keyLandmark_list.get(t1_idx).get(pt).getL() ) );
				}
				else if(t2 == t1){
					for(int pt=0; pt<npt; pt++)
						point_list.add( new Point( keyLandmark_list.get(t1_idx).get(pt).getL() ) );
		 		}
				else{
					for(int pt=0; pt<npt; pt++){ 
						p1[pt] = keyLandmark_list.get(t1_idx).get(pt).getL(); // 1st indices for point, 2nd is for coordinate
						p2[pt] = keyLandmark_list.get(t2_idx).get(pt).getL();
						double[] p = new double[nd];
						for(int d=0; d<nd; d++)
							p[d] = p1[pt][d] + (p2[pt][d]-p1[pt][d])*(t-t1)/(t2-t1);
						point_list.add( new Point( p ) );
						//if (pt==0)
						//	IJ.log("x = "+ point_list.get(0).getL()[0] + " ; y = "+ point_list.get(0).getL()[1] );
					}
				}
				Landmark_list.add( point_list );
				//IJ.log("x = "+ Landmark_list.get(t-1).get(0).getL()[0] + " ; y = "+ Landmark_list.get(t-1).get(0).getL()[1] );	
				//IJ.log("list size "+Landmark_list.size());
			}

			// display successively the interpolated line on the image (for debug)
			//float[] pp1=new float[nd], pp2=new float[nd];
			/*for(int t=1; t<=nFrames; t++)
			{	
				float[] pp1 = Landmark_list.get(t-1).get(0).getL();
				float[] pp2 = Landmark_list.get(t-1).get(1).getL();
				//IJ.log("x = "+ pp1[0] + " ; y = "+ pp1[1]);
				Line line = new Line((double)pp1[0],(double)pp1[1],(double)pp2[0],(double)pp2[1]);
				imp.getProcessor().draw(line);
				imp.updateAndDraw();
				IJ.log("line length at " + t +" : " + line.getLength());
			}
			*/

			// get the current time point as the reference frame
			int ref_t =(int)(imp.getCurrentSlice()/slice_per_frame)+1;
			IJ.log("reference frame : " + ref_t);
			ArrayList<Point> ref_points = Landmark_list.get(ref_t-1);
			
			ImageStack stack = new ImageStack(inputDims[0],inputDims[1]);

			
			int nCZT = inputDims[2]*inputDims[3]*inputDims[4];
			IJ.log("reference slice :"+ref_t);
			
			IJ .log("Processing drift corrected sequence ...");
			
			for(int idx=1; idx<=nCZT; idx++)
			{
			//int idx = 179;
				//SimilarityModel2D model = new SimilarityModel2D();
				TranslationModel2D model = new TranslationModel2D();
				InverseCoordinateTransform ict = model;
				Mapping< InverseCoordinateTransform > mapping = new InverseTransformMapping< InverseCoordinateTransform >( ict );
			
				int t = (int)Math.ceil((double)idx/slice_per_frame);
				//IJ.log("warping frame "+t);
				ImageProcessor ipSource = imp.getStack().getProcessor(idx);
				ImageProcessor ipTarget = ipSource.createProcessor( inputDims[0],inputDims[1] );
				ArrayList< PointMatch > matches = new ArrayList< PointMatch >();
				for ( int pt = 0; pt < npt; ++pt )
				{
					matches.add( new PointMatch( Landmark_list.get(t-1).get(pt) , ref_points.get(pt) ) );
					//IJ.log("pt : "+ Landmark_list.get(t-1).get(pt).getL()[0] + " , "+ref_points.get(pt).getL()[0]);
					//IJ.log("pt : "+ Landmark_list.get(t-1).get(pt).getL()[1] + " , "+ref_points.get(pt).getL()[1]);
					
				}
				try
				{	model.fit( matches );
				}
				catch ( final NotEnoughDataPointsException e )
				{	IJ.showMessage( "Not enough landmarks selected to find a transformation model." );
					return;
				}
				//catch ( final IllDefinedDataPointsException e )
				//{	IJ.showMessage( "The set of landmarks is ill-defined in terms of the desired transformation." );
				//	return;
				//}
				if ( showMatrix ){
					final double[] flatmatrix = new double[6];
					( ( Affine2D< ? > )model ).toArray( flatmatrix );
					IJ.log("Matrix: " + Arrays.toString(flatmatrix));
				}
				if ( interpolate ){
					ipSource.setInterpolationMethod( ImageProcessor.BILINEAR );
					mapping.mapInterpolated( ipSource, ipTarget );
				}
				else{
					mapping.map( ipSource, ipTarget );
				}
				stack.addSlice(ipTarget);				
			}
			
			ImagePlus imp_transform = new ImagePlus("transform sequence", stack);
			// imp_transform.show();
			
			imp_transform.setDimensions(inputDims[2], inputDims[3], inputDims[4]);
			imp_transform.setOpenAsHyperStack(true);
			imp_transform.setCalibration(imp.getCalibration());
			CompositeImage comp_imp_transform = new CompositeImage(imp_transform,CompositeImage.COMPOSITE);
			comp_imp_transform.show();
			
			IJ.log("-------------Drift Correction Done! --------------");	

			// display the new sequence
		}
		
		public static void main(final String... args) {
			
			new ij.ImageJ();
			IJ.open( "A_sequence)_file.tiff" );
			IJ.getImage().setT(100);
			//IJ.run(imp, "Z Project...", "projection=[Max Intensity] all");
			new Manual_Registration().run("");
		}
	}
