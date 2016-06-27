import java.util.ArrayList;


public class ppsdreal{
	protected int [][] kMap;//mark of detected synapses from 1 to nSyn0
	protected double [][] zMap;//zscores of detected synapses
	protected double [][] thrMap;//threshold map (10,20,...250)
	protected double thrZ;//lowest zscore of detected synapses
	protected int nSyn0;//number of detected synapses
	
    public ppsdreal(double [][] G,double [][] Gt, boolean[][] kMask, ParaP p, ParaQ q) {
    	
    	kMap = new int[Gt.length][Gt[0].length];
    	zMap = new double[Gt.length][Gt[0].length];
    	thrMap = new double[Gt.length][Gt[0].length];
    	//initialization
    	boolean doneAll = false;
    	int loopCnt = 0;
    	nSyn0 = 1;
    	//double var0 = q.var;
    	BasicMath BM = new BasicMath();

    	scanAll3 scanMap = new scanAll3(G,Gt,kMask,p,q);// thresholding images with multi thresholds(50:10:220)
    	ArrayList<Double> Iter_PopUplist = new ArrayList<Double>(); //zscores of detected synapses
    	//
    	while(!doneAll){
    		loopCnt = loopCnt + 1;
    		ArrayList<Integer[]> idxUpdt = new ArrayList<Integer[]>();
    		BestSyn bestSynMap = new BestSyn(scanMap, p, Iter_PopUplist);//find best region with highest zscore
    		Iter_PopUplist = new ArrayList<Double>(bestSynMap.PopUplist);
    		// further scan within best region -----
    		if(BM.Allfalse(bestSynMap.kMap0)){
    			doneAll = true;
    		}
    		else{//start scan, find if there is a region score higher inside the best region
    			SynInBest SynBBest = new SynInBest(Gt,bestSynMap,p,q);
    			if(SynBBest.foundOne){
    				//insert the found region(marked in SynBBest) into kMap, zMap, thrMap
    				for(int i=0;i<SynBBest.kMap2.length;i++){
    					for(int j=0;j<SynBBest.kMap2[0].length;j++){
    						if(SynBBest.kMap2[i][j]){
		    		    	            kMap[i][j] = nSyn0;
		    		    	            zMap[i][j] = SynBBest.zMap2[i][j];
		    		    	            thrMap[i][j] = SynBBest.thrMap2[i][j];
    						}
    					}
    				}
    			}
    			else{//otherwise, insert the original region(marked in bestSynMap) into kMap, zMap, thrMap
    				for(int i=0;i<bestSynMap.kMap0.length;i++){
    					for(int j=0;j<bestSynMap.kMap0[0].length;j++){
    						if(bestSynMap.kMap0[i][j]){
		    		    	            kMap[i][j] = nSyn0;
		    		    	            zMap[i][j] = bestSynMap.zMap0[i][j];
		    		    	            thrMap[i][j] = bestSynMap.thrMap0[i][j];
    						}
    					}
    				}
    			}
    			nSyn0++;
    			//update kMask: kMask(kMap>0) = 0;
    			//update idxUpdt: idxUpdt = find(kMap2>0);
    			for(int i=0;i<kMask.length;i++){
    				for(int j=0;j<kMask[0].length;j++){
    					if(kMap[i][j]>0)
    						kMask[i][j] = false;
    					if(SynBBest.kMap2[i][j])
    						idxUpdt.add(new Integer[] {i,j});
    				}
    			}
    			scanMap.scanUpdtCrop(G,Gt,kMask,idxUpdt,p,q);
    		}
    	}
    	nSyn0--;
    }
    // post processing with size constrain and pvalue constrain
    public boolean [][] ppsd_post(double [][] Gt, ParaP p, ParaQ q){
    	System.out.println("Analysis----");
    	ImageHandling IH = new ImageHandling();

    	boolean[][] kMapi = new boolean[kMap.length][kMap[0].length];
    	
    	double thrZx = thrZ<-2? -2:thrZ; //pvalue threshold should not be too small, otherwise FDR control will be not useful
    	int[] roi_size = new int[nSyn0];
    	for(int i=0;i<kMap.length;i++)
    		for(int j=0;j<kMap[0].length;j++)
    			if(kMap[i][j]!=0)
    				roi_size[kMap[i][j]-1]++;
    	
    	for(int jj=0;jj<nSyn0;jj++){
    		int idxSyn0 = jj+1;
    		int[][] idx = new int[roi_size[jj]][2];
    		int idx_cnt = 0;
        	for(int i=0;i<kMap.length;i++)
        		for(int j=0;j<kMap[0].length;j++)
        			if(kMap[i][j]==idxSyn0){
        				idx[idx_cnt][0] = i;
        				idx[idx_cnt++][1] = j;
        			}

        	if(idx.length<10*q.Pix_per_synapse){
        		boolean[][] bmask  = new boolean[q.Ny][q.Nx];
        		boolean[][] tmask  = new boolean[q.Ny][q.Nx];
				for (int i = 0; i < q.Ny; i++) {
					for (int j = 0; j < q.Nx; j++) {
						bmask[i][j] = true;
						tmask[i][j] = false;
					}
				}
				double zScore = IH.scanOneSyn(Gt,tmask,bmask,idx,p,q);
        		if(zScore>thrZx){
        			for(int i=0;i<idx.length;i++){
        				kMapi[idx[i][0]][idx[i][1]] = true;
        			}
        			
        		}
        	}
    	}
    	
    	return kMapi;
    }
   

}
