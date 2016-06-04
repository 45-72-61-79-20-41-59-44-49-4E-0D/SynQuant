import java.awt.Color;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

import javax.imageio.ImageIO;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import ij.*;
import ij.gui.*;
import ij.measure.ResultsTable;
import ij.plugin.*;
import ij.process.*;

/**
 * This plugin 
 * 1. segments synapse/puncta using Probability principled synapse detection
 * algorithm modified from "PPSD: Probability Principled Synapse Detection"
 * 2. extract dendrite from dendrite channel and use linear regression to anaylyse
 * relationships between synapse and 
 * 
 * @author Congchao Wang
 * @contact ccwang@vt.edu
 * @version 1.0
 * @date 2016-05-31
 *
 *
 */

public class SynQuant_ implements PlugIn {
	// image data
	protected ImageStack stack; // Current ImagePlus stack
	protected ImagePlus imp; //Current Image
	protected double[][] G; //2D matrix saving gray scale Image data
	protected double[][] Gt; //stabilized image
	protected int ntry=1; //decide max number of neighbors considered in zscore calculation
	protected int width; //image width
	protected int height; //image height
	protected boolean[][] kSynR1; //bianry synapse Map
	protected int syn_chl;
	protected int den_chl;
	protected double fdr;

	public void run(String arg) {
		try {
			if (showDialog())
				synQuant_real();
		} catch (IndexOutOfBoundsException e) {
			e.printStackTrace();
		}
	}
	public boolean showDialog() 
	{
		// Get input parameter
		GenericDialog gd = new GenericDialog("SynQuant - Data and Parameter Setting");
		gd.addNumericField("Synapse Channel (if NA, input 0)#: ", 2, 1);//2.5-3.5 good
		gd.addNumericField("Dendrite Channel (if NA, input 0)#: ", 3, 1);//2.5-3.5 good
		gd.addNumericField("FDR Control: ", 0.05, 3);//2.5-3.5 good
		gd.showDialog();
		if (gd.wasCanceled()){
			return false;
		}
		syn_chl = (int) gd.getNextNumber();
		den_chl = (int) gd.getNextNumber();
		fdr = gd.getNextNumber();
		if(Double.isNaN(fdr) || (fdr<=0) || (fdr>=1)){
			IJ.showMessage("Invalid parameter(s).\n" + "0-1 needed.");
			return false;
		}
		return true;
	}
	public void synQuant_real() {
		imp = WindowManager.getCurrentImage();
		stack = imp.getStack();
		width = imp.getWidth();
		height = imp.getHeight();
		
		double vox_x = imp.getCalibration().pixelWidth; 
		if(vox_x==1)//simulated data
			vox_x = 2.0757e-7;
		else
			vox_x = vox_x*1e-6;//real data
		
		int type = imp.getType();
		if(syn_chl>0 & den_chl>0){
			//synapse detection
			short[] synArr = stack2array(type, stack, syn_chl);
			ppsd syn_det = new ppsd(synArr, width, height,vox_x, fdr);
			//dendrite extraction
			short[] denArr = stack2array(type, stack, den_chl);
			GrowNeurite den_det = new GrowNeurite(denArr, width, height,vox_x);
			
			LinearTest LT = new LinearTest(syn_det,den_det);
			//show features
			ResultsTable Ft_table = new ResultsTable();
			int DenCnt = 0;
			for (int i=0;i<den_det.CurNum;i++) {
				Ft_table.incrementCounter();
				DenCnt = DenCnt+1;
				Ft_table.addLabel("Dendrite piece #"+DenCnt);
				Ft_table.addValue("Number of synapse on it",LT.denSynNum[i]);
				Ft_table.addValue("Dendrite length",LT.DenFeatures[i][0]);
				Ft_table.addValue("Dendrite scale",LT.DenFeatures[i][1]);
				Ft_table.addValue("Dendrite intensity",LT.DenFeatures[i][2]);
			}
			Ft_table.showRowNumbers(false);
			Ft_table.show("Feature Table");
			
			ResultsTable sm_table = new ResultsTable();

			sm_table.incrementCounter();
			sm_table.addLabel("Coefficents:");
			sm_table.addValue("Length",LT.betas[0]);
			sm_table.addValue("Scale",LT.betas[1]);
			sm_table.addValue("Intensity",LT.betas[2]);
			sm_table.addValue("Intercept",LT.betas[3]);
			/*add pvalues of regression results
			 * sm_table.incrementCounter();
			sm_table.addLabel("P-values:");
			sm_table.addValue("Length",LT.pvalues[0]);
			sm_table.addValue("Scale",LT.pvalues[1]);
			sm_table.addValue("Intensity",LT.pvalues[2]);
			sm_table.addValue("Intercept",LT.pvalues[3]);*/
			sm_table.showRowNumbers(false);
			sm_table.show("Summary table");
			
		}
		else{
			if(syn_chl>0){
				short[] synArr = stack2array(type, stack, syn_chl);
				ppsd syn_det = new ppsd(synArr, width, height,vox_x, fdr);
			}
			if(den_chl>0){
				short[] denArr = stack2array(type, stack, den_chl);
				GrowNeurite den_det = new GrowNeurite(denArr, width, height,vox_x);
			}
		}

	}
	public short[] stack2array(int type, ImageStack stack, int channel){
		int mask=0xff;
		int nPixels=width*height;
		short[] imArray = new short[nPixels];
		if (type == ImagePlus.GRAY16)
		{
			mask=0xffff;
			//Find min, max. Copy to imArray
			short[] pixels = (short[])stack.getPixels(channel);
			int intP = (int)(mask&pixels[0]);
			for (int i=0; i<nPixels; i++)
			{
				intP=(int)(mask&pixels[i]);
				short p = (short)(intP/2);
				imArray[i]=p;
			}
		}  
		else if (type == ImagePlus.GRAY8) 
		{
			mask=0xff;
			//Find min, max. Copy to imArray
			byte[] pixels = (byte[])stack.getPixels(channel);
			for (int i=0; i<nPixels; i++)
			{
				short p=(short)(mask&pixels[i]);
				imArray[i]=p;
			}
		}
		else
		{
			IJ.log("Pixel format not supported");
			return null;
		}
		return imArray;
	}
	/**
	 * Main method for debugging.
	 *
	 * For debugging, it is convenient to have a method that starts ImageJ,
	 * loads an image and calls the plugin, e.g. after setting breakpoints.
	 *
	 * @param args
	 *            unused
	 */
	/*public static void main(String[] args) {
		// set the plugins.dir property to make the plugin appear in the Plugins
		// menu
		Class<?> clazz = SynQuant_.class;
		String url = clazz.getResource("/" + clazz.getName().replace('.', '/') + ".class").toString();
		String pluginsDir = url.substring(5, url.length() - clazz.getName().length() - 6);
		System.setProperty("plugins.dir", pluginsDir);

		// start ImageJ
		new ImageJ();

		// open the Clown sample
		//ImagePlus image = IJ
		//		.openImage("C:\\Users\\ccwang\\Desktop\\DS2_3_1024_0005_001_0.tiff");//C2-3-weak-z_Maximum intensity projection.tif");
		//image.show();

		// run the plugin
		//IJ.runPlugIn(clazz.getName(), "");
	}*/

	void error() {
		IJ.showMessage("PPSD", "Error");
	}
}
