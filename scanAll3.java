import java.util.ArrayList;
import java.util.Collections;

public class scanAll3{
	protected int [][][] kMapx;//regions numbering for each threshold
	protected double [][][] zMapx;//zscore maps for each threshold
	protected int Nthr; //number of thresholds used
	protected int[] thrrg; // thresholds

	//SCANALL scan all thresholds and build the zscore map of each threshold
	//Parameters are initialization for the first iteration
	//kMask is binary synapse map with background 1 and detected synaspe 0
	public scanAll3(double[][] G, double[][] Gt, boolean[][] kMask, ParaP p, ParaQ q) {
		
		Nthr = (p.thr1-p.thr0)/p.thrg +1;
		thrrg = new int[Nthr];
		for(int i = 0;i<Nthr;i++)
			thrrg[i] = p.thr0+i*p.thrg; 
		ImageHandling IH = new ImageHandling();
		kMapx = new int[Nthr][][];
		zMapx = new double[Nthr][][];
		for(int ii=0;ii<Nthr;ii++){
			int thr = thrrg[ii];
			boolean [][] K1 = thresholding(Gt, thr,kMask);

			boolean [][] K1b = IH.imMorpho(K1, kMask, p.min_size);
			int [][] K1b_cc = IH.bwlabel(K1b,8);
			int ccN = IH.NextLabel;				
			kMapx[ii] = new int[Gt.length][Gt[0].length];
			zMapx[ii] = new double[Gt.length][Gt[0].length];
			if(ccN==0)
			{
				continue;
			}
			double [] reg_mean = new double[ccN];
			int [] reg_size = new int[ccN];
			double [] z_score = new double[ccN];
			//int [] k_x = new int[ccN]; // region count
			int [][] reg_cor = new int[ccN][4]; //ymin xmin ymax xmax  
			for(int i=0;i<ccN;i++){
				reg_cor[i][0] = G.length;
				reg_cor[i][1] = G[0].length;
			}
			
			for(int i=0;i<K1b_cc.length;i++){
				for(int j=0;j<K1b_cc[0].length;j++){
					if(K1b_cc[i][j]!=0){
						reg_mean[K1b_cc[i][j]-1] += G[i][j];
						reg_size[K1b_cc[i][j]-1]++;
						if(i<reg_cor[K1b_cc[i][j]-1][0])
							reg_cor[K1b_cc[i][j]-1][0]=i;
						if(j<reg_cor[K1b_cc[i][j]-1][1])
							reg_cor[K1b_cc[i][j]-1][1]=j;
						if(i>reg_cor[K1b_cc[i][j]-1][2])
							reg_cor[K1b_cc[i][j]-1][2]=i;
						if(j>reg_cor[K1b_cc[i][j]-1][3])
							reg_cor[K1b_cc[i][j]-1][3]=j;
					}
				}
			}

			for(int i=0;i<ccN;i++){
				int[][] reg0 = new int[reg_size[i]][2];
				int reg0_cnt = 0;
				for(int y=0;y<K1b_cc.length;y++){
					for(int x=0;x<K1b_cc[0].length;x++){
						if(K1b_cc[y][x]==i+1){
							reg0[reg0_cnt][0] = y;
							reg0[reg0_cnt++][1] = x;
						}
				}
				}
				if(reg_size[i]>p.max_size || reg_size[i]<p.min_size){
					z_score[i] = 0;
					//continue;
				}
				else if(reg_mean[i]/reg_size[i]<p.minIntensity){
					z_score[i] = 0;
					//continue;
				}
				else{
					double LH = reg_cor[i][3]-reg_cor[i][1]+1;
					double LW = reg_cor[i][2]-reg_cor[i][0]+1;
					double ratio = LH>LW? LH/LW: LW/LH;
					if(ratio>p.maxWHratio || (double)reg_size[i]/(LH*LW)<p.minfill){
						z_score[i] = 0;
						//continue;
					}
					else
						z_score[i] = IH.scanOneSyn(Gt,K1,kMask,reg0,p,q);
						System.out.println("Threshold "+ii+"-region:"+i+"-zscore:"+z_score[i]);
				}
				for(int j=0;j<reg0.length;j++){
					kMapx[ii][reg0[j][0]][reg0[j][1]] = i+1;
				}
				if(z_score[i]!=0){
					for(int j=0;j<reg0.length;j++)
						zMapx[ii][reg0[j][0]][reg0[j][1]] = z_score[i];
				}
			}
		}
	}
	//SCANALL scan all thresholds and find the best map in a cropped region
	public scanAll3/*scanAll3Crop*/(double[][] Gt, boolean[][] Kmask,ParaP p,ParaQ q, int thr0) {//find the best synapse region
		
		Nthr = (p.thr1-p.thr0)/p.thrg +1;
		thrrg = new int[Nthr];
		for(int i = 0;i<Nthr;i++)
			thrrg[i] = p.thr0+i*p.thrg; 
		ImageHandling IH = new ImageHandling();
		kMapx = new int[Nthr][][];
		zMapx = new double[Nthr][][];
		int Ny = Gt.length;
		int Nx = Gt[0].length;
		int ofst = 5;
		for(int ii=0;ii<Nthr;ii++){
			int thr = thrrg[ii];
			kMapx[ii] = new int[Ny][Nx];
			zMapx[ii] = new double[Ny][Nx];
			int [][] Kx = kMapx[ii];
			double [][] Zx = zMapx[ii];
			if(!(thr>thr0)){
				continue;
			}
			boolean [][] K = thresholding(Gt, thr,Kmask);
			//crop
			ArrayList<Integer> r = new ArrayList<Integer>();
			ArrayList<Integer> c = new ArrayList<Integer>();
			for(int i=0;i<K.length;i++){
				for(int j=0;j<K[0].length;j++){
					if(K[i][j]){
						r.add(i);
						c.add(j);
					}
				}
			}

			if(r.size()==0 || c.size()==0){
				continue;
			}
				
			//Collections.sort(c);
			int rgr_min = Math.max(0, Collections.min(r)-ofst);
			int rgr_max = Math.min(K.length-1, Collections.max(r)+ofst);
			int rgc_min = Math.max(0, Collections.min(c)-5);
			int rgc_max = Math.min(K[0].length-1, Collections.max(c)+5);
			double[][] Gt1 = new double[rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			boolean[][] Kmask1 = new boolean[rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			boolean[][] K1 = new boolean[rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			double[][] Zx1 = new double[rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			int[][] Kx1 = new int[rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			
			for(int i=0;i<Gt1.length;i++){
				for(int j=0;j<Gt1[0].length;j++){
					Gt1[i][j] = Gt[rgr_min+i][rgc_min+j];
					Kmask1[i][j] = Kmask[rgr_min+i][rgc_min+j];
					K1[i][j] = K[rgr_min+i][rgc_min+j];
					Zx1[i][j] = Zx[rgr_min+i][rgc_min+j];
					Kx1[i][j] = Kx[rgr_min+i][rgc_min+j];
				}
			}
	
			boolean [][] K1b = IH.imMorpho(K1, Kmask1, p.min_size);
			int [][] K1b_cc = IH.bwlabel(K1b,8);
			int ccN1 = IH.NextLabel;
			if(ccN1==0)
			{
				continue;
			}
			int [] reg_size = new int[ccN1];
			double [] z_score = new double[ccN1];

			
			for(int i=0;i<K1b_cc.length;i++){
				for(int j=0;j<K1b_cc[0].length;j++){
					if(K1b_cc[i][j]!=0){
						reg_size[K1b_cc[i][j]-1]++;
					}
				}
			}
			for(int i=0;i<ccN1;i++){
				int[][] smaskIdx = new int[reg_size[i]][2];
				int smaskIdx_cnt = 0;
				for(int y=0;y<K1b_cc.length;y++){
					for(int x=0;x<K1b_cc[0].length;x++){
						if(K1b_cc[y][x]==i+1){
							smaskIdx[smaskIdx_cnt][0] = y;
							smaskIdx[smaskIdx_cnt++][1] = x;
						}
				}
				}
				
				if(reg_size[i]>p.max_size || reg_size[i]<p.min_size){
					z_score[i] = 0;
					//continue;
				}
				else{
					z_score[i] = IH.scanOneSyn(Gt1,K1b,Kmask1,smaskIdx,p,q);
				}
				for(int j=0;j<smaskIdx.length;j++){
					kMapx[ii][smaskIdx[j][0]+rgr_min][smaskIdx[j][1]+rgc_min] = j+1;
				}
				if(z_score[i]!=0)
					for(int j=0;j<smaskIdx.length;j++)
						zMapx[ii][smaskIdx[j][0]+rgr_min][smaskIdx[j][1]+rgc_min] = z_score[i];
			}
		}
	}
	//update the synapse map and zscore map
	//bmask is the background mask
	//idxUpdt is the pixel positions of detected synapse in form iteration
	public void scanUpdtCrop(double[][] G, double[][] Gt, boolean[][] bmask, ArrayList<Integer[]> idxUpdt, ParaP p,
			ParaQ q) {
		
		//q.ntry = 1 here, but what is the right range of neighboring size?
		ImageHandling IH = new ImageHandling();
		BasicMath mBM = new BasicMath();
		for(int ii=0;ii<Nthr;ii++){
			int[][] Kx = kMapx[ii];
			double [][] Zx = zMapx[ii];
			int[] u0 = new int [idxUpdt.size()];
			for(int i=0;i<idxUpdt.size();i++){
				u0[i] = Kx[idxUpdt.get(i)[0]][idxUpdt.get(i)[1]];
			}
			int lu0 = mBM.unique(u0).length;
			//System.out.println("number of idx "+lu0);
			if(lu0==1 && mBM.vectorMax(u0)==0){
				continue;
			}
			if(lu0>1){
				for(int i=0;i<idxUpdt.size();i++){
					Kx[idxUpdt.get(i)[0]][idxUpdt.get(i)[1]]=0;
					Zx[idxUpdt.get(i)[0]][idxUpdt.get(i)[1]]=0;
				}
				continue;
			}
			int idx0 = (int) Math.round((double)mBM.vectorSum(u0)/u0.length);//there should be only one value in u0
			boolean [][] maskUpdt = new boolean[Kx.length][Kx[0].length];
			for(int i=0;i<Kx.length;i++){
				for(int j=0;j<Kx[0].length;j++){
					if(Kx[i][j]==idx0){
						maskUpdt[i][j] = true;
						Kx[i][j] = 0;//
						Zx[i][j] = 0;//
					}
					else
						maskUpdt[i][j] = false;
				}
			}
			int nSyn0 = mBM.matrix2DMax(Kx);
			int thr = thrrg[ii];
			boolean [][] K = thresholding(Gt, thr,maskUpdt);
			//crop
			ArrayList<Integer> r = new ArrayList<Integer>();
			ArrayList<Integer> c = new ArrayList<Integer>();
			for(int i=0;i<K.length;i++){
				for(int j=0;j<K[0].length;j++){
					if(K[i][j]){
						r.add(i);
						c.add(j);
					}
				}
			}
			if(r.size()==0 || c.size()==0){
				continue;
			}
			
			int ofst = 5;// cut the region with width+10 * height+10
			int rgr_min = Math.max(0, Collections.min(r)-ofst);
			int rgr_max = Math.min(K.length-1, Collections.max(r)+ofst);
			int rgc_min = Math.max(0, Collections.min(c)-5);
			int rgc_max = Math.min(K[0].length-1, Collections.max(c)+5);
			double[][] Gt1 = new double[rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			boolean[][] bmask1 = new boolean[rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			boolean[][] K1 = new boolean[rgr_max-rgr_min+1][rgc_max-rgc_min+1];
			
			for(int i=0;i<Gt1.length;i++){
				for(int j=0;j<Gt1[0].length;j++){
					Gt1[i][j] = Gt[rgr_min+i][rgc_min+j];
					bmask1[i][j] = bmask[rgr_min+i][rgc_min+j];
					K1[i][j] = K[rgr_min+i][rgc_min+j];
				}
			}

			boolean [][] tmask1 = IH.imMorpho(K1, bmask1, p.min_size);
			int [][] tmask1_cc = IH.bwlabel(tmask1,8);
			int ccN1 = IH.NextLabel;
			if(ccN1==0)
			{
				continue;
			}
			double [] reg_mean = new double[ccN1];
			int [] reg_size = new int[ccN1];
			double [] z_score = new double[ccN1];
			//int [] k_x = new int[ccN]; // region count
			int [][] reg_cor = new int[ccN1][4]; //ymin xmin ymax xmax  
			for(int i=0;i<ccN1;i++){
				reg_cor[i][0] = tmask1_cc.length;
				reg_cor[i][1] = tmask1_cc[0].length;
			}
			
			for(int i=0;i<tmask1_cc.length;i++){
				for(int j=0;j<tmask1_cc[0].length;j++){
					if(tmask1_cc[i][j]!=0){
						reg_mean[tmask1_cc[i][j]-1] += G[i][j];
						reg_size[tmask1_cc[i][j]-1]++;
						if(i<reg_cor[tmask1_cc[i][j]-1][0])
							reg_cor[tmask1_cc[i][j]-1][0]=i;
						if(j<reg_cor[tmask1_cc[i][j]-1][1])
							reg_cor[tmask1_cc[i][j]-1][1]=j;
						if(i>reg_cor[tmask1_cc[i][j]-1][2])
							reg_cor[tmask1_cc[i][j]-1][2]=i;
						if(j>reg_cor[tmask1_cc[i][j]-1][3])
							reg_cor[tmask1_cc[i][j]-1][3]=j;
					}
				}
			}
			for(int i=0;i<ccN1;i++){
				int[][] smaskIdx = new int[reg_size[i]][2];
				int smaskIdx_cnt = 0;
				for(int y=0;y<tmask1_cc.length;y++){
					for(int x=0;x<tmask1_cc[0].length;x++){
						if(tmask1_cc[y][x]==i+1){
							smaskIdx[smaskIdx_cnt][0] = y;
							smaskIdx[smaskIdx_cnt++][1] = x;
						}
				}
				}
				
				if(reg_size[i]>p.max_size || reg_size[i]<p.min_size){
					z_score[i] = 0;
					//continue;
				}
				else if(reg_mean[i]/reg_size[i]<p.minIntensity){
					z_score[i] = 0;
					//continue;
				}
				else{
					int LH = reg_cor[i][3]-reg_cor[i][1]+1;
					int LW = reg_cor[i][2]-reg_cor[i][0]+1;
					double ratio = LH>LW? LH/LW: LW/LH;
					if(ratio>p.maxWHratio || reg_size[i]/(LH*LW)<p.minfill){
						z_score[i] = 0;
						//continue;
					}
					else{
						z_score[i] = IH.scanOneSyn(Gt1,tmask1,bmask1,smaskIdx,p,q);
					}
				}
				for(int j=0;j<smaskIdx.length;j++){
					kMapx[ii][smaskIdx[j][0]+rgr_min][smaskIdx[j][1]+rgc_min] = j+1+nSyn0;
				}
				if(z_score[i]!=0)
					for(int j=0;j<smaskIdx.length;j++)
						zMapx[ii][smaskIdx[j][0]+rgr_min][smaskIdx[j][1]+rgc_min] = z_score[i];
				System.out.println(i+"Updated z_score："+z_score[i]);
			}
		}
	}
//threshold the image with a threshold and mask
	public boolean [][] thresholding(double[][] Gt, int thr, boolean [][] kMask){
		boolean [][] K1 = new boolean [Gt.length][Gt[0].length];
		for(int i=0;i<Gt.length;i++)
			for(int j=0;j<Gt[0].length;j++)
				K1[i][j] = (Gt[i][j]>thr) && kMask[i][j];
		return K1;
	}
}