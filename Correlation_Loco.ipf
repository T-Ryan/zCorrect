#pragma rtGlobals=1		// Use modern global access method and strict wave access.
Menu "Rozan_Procedures"
	"Correlation and Significance test for a time-series |e.g. locomotion read-out| vs Calcium traces", /Q,  CorrSign()
End
End





//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
//Modified TR-01/05/19 modified the way the distribution of random corr values is used to compute pvalue (previously only passed negative correlations |||||||||||||||||||||||||||||||//
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||//
Function CorrSign2()

	////// GUI ///////
	string MatrixWaveName,LocoWaveName,MWN,LWN
	string SPEARorPEAR_select,histograph
	variable Hz=10
	LocoWaveName=wavelist("*", ";", "DIMS:1")
	MatrixWaveName=wavelist("*", ";", "DIMS:2")
	prompt SPEARorPEAR_select, "SPEAR or PEAR?", popup "Spearman's Rank Correlation;Pearson's Correlation (assumes normally distributed data)"
	prompt histograph, "Would you like me to plot a histogram of the correlation coefficient values?", popup "Yes please :);No, I'm good. Thanks though :)"
	prompt LWN, "Select a 1D wave time-series", popup LocoWaveName
	prompt MWN, "Select a 2D wave containing the calcium traces of your ROIs", popup MatrixWaveName
	prompt Hz, "What is the acquisition frequency?"
	doprompt "Correlation of time-series + significance test and mult. comp. corr.",SPEARorPEAR_select,histograph,LWN,MWN,Hz
	if(V_flag==1)
		abort
	endif

	wave stimcycle
	variable sigpvalue = 0.05
	variable nshuffles = 1000
	variable pvalue = 0


	print SPEARorPEAR_select
	variable SPEARorPEAR
	if (stringmatch(SPEARorPEAR_select,"Spearman's Rank Correlation"))
		SPEARorPEAR=0
	elseif (stringmatch(SPEARorPEAR_select,"Pearson's Correlation (assumes normally distributed data)"))
		SPEARorPEAR=1
	endif
	
	duplicate/FREE/o $MWN All
	duplicate/FREE/o $LWN Loco

	wave stim
	wave loco_timex

	// make a new folder "Correlation" if the folder does not aready exists
	string Correlation
	Correlation="root:Correlation"
	if (!datafolderexists(Correlation))
		newdatafolder root:Correlation
	elseif (datafolderexists(Correlation))
		print "The folder 'Correlation' already exists"
	endif
	variable NumFrames=dimsize(Loco,0)
	make/FREE/o/n=(dimsize(All,1)) List_Corr_LOCO_CCvalues
	make/FREE/o/n=(dimsize(All,1)) List_Corr_LOCO_Pvalues

	Make/FREE/o/n=(NumFrames) RandomWave 
	Make/FREE/o/n=(NumFrames) RandomCorValue

	Make/FREE/o/n=0 ROI_list_Corr_NS
	Make/FREE/o/n=0 ROI_list_Corr_pos
	Make/FREE/o/n=0 ROI_list_Corr_neg

	variable LenRandom=dimsize(RandomWave,0)
	variable LenAll=dimsize(All,0)
	variable LenROIs=dimsize(All,1)

	// Generating a wave with a 1000 random starting points for the circular shifting
	RandomWave=0
	variable i,j,k,ii,jj,kk,iii,jjj,count1
	variable willekeur
	for (j=0;j<LenRandom;j+=1)
		i=0
		do
			willekeur=round(enoise(LenAll/2))+(LenAll/2)
			Findvalue/V=(willekeur) RandomWave
			If (V_Value==-1 && willekeur!=0)
				RandomWave[j]=willekeur
				i=1
			endif
		while (i!=1)
	endfor

	// Big Loop through all ROIs

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	// SPEARMAN                                                                                                      //
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	if (SPEARorPEAR==0)

		variable FirstValue
		variable HighValue
		//Correlation: Loop through all ROIs
		FOR (i=0;i<LenROIs;i+=1)
			count1=0
			variable CorrList
			duplicate/FREE/o/r=[][i] All temp1
			StatsRankCorrelationTest/ALPH=0.05/T=1/Q Loco, temp1
			wave W_StatsRankCorrelationTest
			CorrList=W_StatsRankCorrelationTest[4]
			List_Corr_LOCO_CCvalues[i]=CorrList
	
			//Correlation done with random shifts a number of times equal to RandomWave size
			for (j=0;j<LenRandom;j+=1)
				duplicate/FREE/o Loco tempStart
				Duplicate/FREE/o/r=[0,RandomWave[j]-1] Loco tempEnd
				Deletepoints/m=0 0,RandomWave[j]+0, tempStart 
				Concatenate/o/np=0 {tempStart,tempEnd}, tempRandom
				if (dimsize(tempRandom,0)!=dimsize(Loco,0)) // Wavelength mismatch
					duplicate/FREE/o Loco tempStart
					Duplicate/FREE/o/r=[0,RandomWave[j]-1] Loco tempEnd
					Deletepoints/m=0 0,RandomWave[j]+1, tempStart // Need to add +1 to make sure there is no wavelength mismatch
					Concatenate/o/np=0 {tempStart,tempEnd}, tempRandom
				endif
				StatsRankCorrelationTest/ALPH=0.05/T=1/Q tempRandom, temp1
				wave W_StatsRankCorrelationTest
				RandomCorValue[j]=W_StatsRankCorrelationTest[4] // 4 for Spearman				This wave now contains correlation values for nshuffle random shuffles
			endfor	
	
			//Calculate the distribution of corr values 
			Make/N=(LenRandom/4)/O Corr_Hist
			Histogram/cum/C/B=1 RandomCorValue, Corr_Hist 
			//Cumulative histogram
			Corr_Hist/=LenRandom										//Divide by number of samples to get CDF	
	
			if(CorrList<0)			
				pvalue = Corr_Hist(List_Corr_LOCO_CCvalues[i])		// Gives pvalue for negative corr values (analagous to the t-value): represents the area under the PDF curve to the LEFT of the corr value 
			else
				pvalue= 1 - Corr_Hist(List_Corr_LOCO_CCvalues[i])	// Take 1-CDF value at corr value (t-value) to acquire the area under the PDF (probability density function) to the RIGHT of the corr-value
			endif
	
			List_Corr_LOCO_Pvalues[i]=pvalue
		ENDFOR

		//stats_BHFDR(List_Corr_LOCO_Pvalues) // Correct for multiple comparison

		// Devide the ROIs in significantly positively or negatively correlated with locomotion, or not significant
		//Changed temporarily to make it more permissive. Original 0.05. Also change in the histogram.
		variable PvalueThr=0.05
		FOR (i=0;i<LenROIs;i+=1)
			if (List_Corr_LOCO_Pvalues[i]>PvalueThr)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_NS,0)),1, ROI_list_Corr_NS
				ROI_list_Corr_NS[dimsize(ROI_list_Corr_NS,0)-1]=i
			elseif (List_Corr_LOCO_Pvalues[i]<=PvalueThr && List_Corr_LOCO_CCvalues[i]>=0)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_pos,0)),1, ROI_list_Corr_pos
				ROI_list_Corr_pos[dimsize(ROI_list_Corr_pos,0)-1]=i
			elseif (List_Corr_LOCO_Pvalues[i]<=PvalueThr && List_Corr_LOCO_CCvalues[i]<0)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_neg,0)),1, ROI_list_Corr_neg
				ROI_list_Corr_neg[dimsize(ROI_list_Corr_neg,0)-1]=i
			endif
		ENDFOR

		string window_name=winname(0,2,1)
		killwindow $window_name
		//killwaves W_StatsRankCorrelationTest, tempRandom, Loco_Corrected

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
		// PEARSON                                                                                                         //
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	elseif (SPEARorPEAR==1)

		FOR (i=0;i<LenROIs;i+=1)
			count1=0
			duplicate/FREE/o/r=[][i] All temp1
			make/FREE/o/n=1 TijdelijkeCorrList
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Determine the highest correlation value by shifting the locomotion time-series minus to pls 2 seconds //
			duplicate/FREE/o Loco tempStart
			Duplicate/FREE/o/r=[(dimsize(Loco,0)-(2*Hz)),dimsize(Loco,0)-1] Loco tempEnd
			Deletepoints/m=0 (dimsize(Loco,0)-(2*Hz)),dimsize(Loco,0), tempStart
			Concatenate/o/np=0 {tempEnd,tempStart}, Loco_Corrected
			StatsLinearCorrelationTest/ALPH=0.05/Q Loco_Corrected, temp1
			wave W_StatsLinearCorrelationTest
			TijdelijkeCorrList[0]=W_StatsLinearCorrelationTest[1]
			for (ii=0;ii<4*Hz;ii+=1)
				FirstValue=Loco_Corrected[0]
				Deletepoints/m=0 0,1, Loco_Corrected
				insertpoints/m=0 dimsize(Loco_Corrected,0),1, Loco_Corrected
				Loco_Corrected[dimsize(Loco_Corrected,0)-1]=FirstValue
				StatsLinearCorrelationTest/ALPH=0.05/Q Loco_Corrected, temp1
				wave W_StatsLinearCorrelationTest
				insertpoints/m=0 dimsize(TijdelijkeCorrList,0),1, TijdelijkeCorrList
				TijdelijkeCorrList[dimsize(TijdelijkeCorrList,0)-1]=W_StatsLinearCorrelationTest[1]
			endfor
			wavestats/q TijdelijkeCorrList
			If (abs(V_max)>abs(V_min))
				HighValue=V_max
			else
				HighValue=V_min
			endif
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			StatsLinearCorrelationTest/ALPH=0.05/Q Loco, temp1
			wave W_StatsLinearCorrelationTest
			List_Corr_LOCO_CCvalues[i]=W_StatsLinearCorrelationTest[1] // 1 for Pearson
			for (j=0;j<nshuffles;j+=1)
				duplicate/FREE/o Loco tempStart
				Duplicate/FREE/o/r=[0,RandomWave[j]-1] Loco tempEnd
				Deletepoints/m=0 0,RandomWave[j], tempStart
				Concatenate/o/np=0 {tempStart,tempEnd}, tempRandom
				if (dimsize(tempRandom,0)!=dimsize(Loco,0)) // Wavelength mismatch
					duplicate/FREE/o Loco tempStart
					Duplicate/FREE/o/r=[0,RandomWave[j]-1] Loco tempEnd
					Deletepoints/m=0 0,RandomWave[j]+1, tempStart // Need to add +1 to make sure there is no wavelength mismatch
					Concatenate/o/np=0 {tempStart,tempEnd}, tempRandom
				endif
				StatsLinearCorrelationTest/ALPH=0.05/T=1/Q tempRandom, temp1
				wave W_StatsLinearCorrelationTest
				RandomCorValue[j]=W_StatsLinearCorrelationTest[1] // 1 for Pearson
			endfor
	
			//Calculate the distribution of corr values 
			Make/N=(LenRandom/4)/O Corr_Hist
			Histogram/CUM/C/B=1 RandomCorValue, Corr_Hist     //Cumulative histogram
			Corr_Hist/=LenRandom
			//Divide by number of samples to get CDF	
			if(CorrList<0)			
				pvalue = Corr_Hist(List_Corr_LOCO_CCvalues[i])		// Gives pvalue for negative corr values (analagous to the t-value): represents the area under the PDF curve to the LEFT of the corr value 
			else
				pvalue= 1 - Corr_Hist(List_Corr_LOCO_CCvalues[i])	// Take 1-CDF value at corr value (t-value) to acquire the area under the PDF (probability density function) to the RIGHT of the corr-value
			endif
	
			List_Corr_LOCO_Pvalues[i]=pvalue
		ENDFOR

		//stats_BHFDR(List_Corr_LOCO_Pvalues) // Correct for multiple comparison

		// Devide the ROIs in significantly positively or negatively correlated with locomotion, or not significant
		FOR (i=0;i<LenROIs;i+=1)
			if (List_Corr_LOCO_Pvalues[i]>0.05)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_NS,0)),1, ROI_list_Corr_NS
				ROI_list_Corr_NS[dimsize(ROI_list_Corr_NS,0)-1]=i
			elseif (List_Corr_LOCO_Pvalues[i]<=0.05 && List_Corr_LOCO_CCvalues[i]>=0)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_pos,0)),1, ROI_list_Corr_pos
				ROI_list_Corr_pos[dimsize(ROI_list_Corr_pos,0)-1]=i
			elseif (List_Corr_LOCO_Pvalues[i]<=0.05 && List_Corr_LOCO_CCvalues[i]<0)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_neg,0)),1, ROI_list_Corr_neg
				ROI_list_Corr_neg[dimsize(ROI_list_Corr_neg,0)-1]=i
			endif
		ENDFOR

		string window_name2=winname(0,2,1)
		killwindow $window_name2
		killwaves W_StatsLinearCorrelationTest, tempRandom, Loco_Corrected

	endif
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Put the new waves in the "Correlation" Folder
	duplicate/o ROI_list_Corr_neg Root:Correlation:ROI_list_Corr_neg
	duplicate/o ROI_list_Corr_pos Root:Correlation:ROI_list_Corr_pos
	duplicate/o ROI_list_Corr_NS Root:Correlation:ROI_list_Corr_NS

	duplicate/o List_Corr_LOCO_CCvalues Root:Correlation:List_Corr_LOCO_CCvalues
	duplicate/o List_Corr_LOCO_Pvalues Root:Correlation:List_Corr_LOCO_Pvalues
	duplicate/o Loco Root:Correlation:Ch3Loco_AvUni
	duplicate/o All Root:Correlation:CorrectedQA_DivMoff

	if (stringmatch(histograph,"Yes please :)"))
		TwoClusterHist()
	endif

END

function stats_BHFDR(w,[ncol])
 
	//-----------------------------------------------------------------
	// parameters
	//-----------------------------------------------------------------
	wave w;
	variable ncol;
 
	//-----------------------------------------------------------------
	// variables
	//-----------------------------------------------------------------
	variable npnts,ndim;
 
	//-----------------------------------------------------------------
	// init variables
	//-----------------------------------------------------------------
 
	if(paramIsDefault(ncol))
		ncol = 0;
	endif
 
	npnts = dimsize(w,0);
	ndim = dimsize(w,1);
 
	//-----------------------------------------------------------------
	// main function
	//-----------------------------------------------------------------
 
	// copy p-values out to new wave and make index wave
	make/free/n=(npnts) wpv,wpi;
 
	if(ndim)
		wpv[] = w[p][ncol];
	else
		wpv[] = w[p];
	endif
 
	// populate index wave
	wpi[] = p;
 
	// sort from highest to lowest
	sort/r wpv,wpv,wpi;
 
	// perform correction
	wpv[] = wpv[p]*(npnts/(npnts-p));
 
	// restore original sorting
	sort wpi,wpv,wpi;
 
	// copy adjusted p-values to original wave
	if(ndim)
		w[][ncol] = wpv[p];
	else
		w[] = wpv[p];
	endif
 
	// set p-values greater than 1 to 1
	w[] = w[p] > 1 ? 1 : w[p];
 
end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This is a function for plotting a histogram for all significant and non-significant correlation coefficient values
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Function TwoClusterHist()

	DFREF saveDFR = GetDataFolderDFR()
	setdatafolder root:Correlation

	wave List_Corr_LOCO_CCvalues
	wave List_Corr_LOCO_Pvalues
	wave ROI_list_Corr_neg
	wave ROI_list_Corr_pos

	Make/FREE/o/n=(dimsize(List_Corr_LOCO_CCvalues,0))  List_CCvalues_All
	Make/o/n=0  List_Sign_CCvalues_All
	variable LenListCC=dimsize(List_Corr_LOCO_CCvalues,0)


	variable i,ii,j,jj,k,kk
	variable PvalueThr=0.05
	For (i=0;i<LenListCC;i+=1)
		List_CCvalues_All[i]=List_Corr_LOCO_CCvalues[i]
		if (List_Corr_LOCO_Pvalues[i]<PvalueThr)
			insertpoints/m=0 dimsize(List_Sign_CCvalues_All,0),1, List_Sign_CCvalues_All
			List_Sign_CCvalues_All[dimsize(List_Sign_CCvalues_All,0)-1]=List_Corr_LOCO_CCvalues[i]
		endif
	Endfor


	////////////// Plot //////////////
	duplicate/FREE/o List_CCvalues_All temp1
	Make/N=80/O Hist_All_ALL;DelayUpdate
	Histogram/B={-1,0.025,80} temp1,Hist_All_ALL

	duplicate/o List_Sign_CCvalues_All temp1
	Make/N=80/O Hist_Sign_All;DelayUpdate
	Histogram/B={-1,0.025,80} temp1,Hist_Sign_ALL

	Display/k=1 Hist_All_ALL
	appendtograph Hist_Sign_ALL
	ModifyGraph mode(Hist_All_ALL)=5,rgb(Hist_All_ALL)=(0,0,0)
	ModifyGraph mode=5,rgb=(0,0,0),hbFill(Hist_Sign_ALL)=2
	Legend/C/N=text1/J/A=RT "\\s(Hist_All_ALL) All correlation coefficient values\r\\s(Hist_Sign_All) Significant correlation coefficient values"

	SetDataFolder saveDFR

END

Function TwoClusterHisttarget(target)

	string target
	string targ= "root:"+target
	DFREF saveDFR = GetDataFolderDFR()
	setdatafolder targ

	wave List_Corr_LOCO_CCvalues
	wave List_Corr_LOCO_Pvalues
	wave ROI_list_Corr_neg
	wave ROI_list_Corr_pos

	Make/FREE/o/n=(dimsize(List_Corr_LOCO_CCvalues,0))  List_CCvalues_All
	Make/o/n=0  List_Sign_CCvalues_All
	variable LenListCC=dimsize(List_Corr_LOCO_CCvalues,0)


	variable i,ii,j,jj,k,kk
	variable PvalueThr=0.05
	For (i=0;i<LenListCC;i+=1)
		List_CCvalues_All[i]=List_Corr_LOCO_CCvalues[i]
		if (List_Corr_LOCO_Pvalues[i]<PvalueThr)
			insertpoints/m=0 dimsize(List_Sign_CCvalues_All,0),1, List_Sign_CCvalues_All
			List_Sign_CCvalues_All[dimsize(List_Sign_CCvalues_All,0)-1]=List_Corr_LOCO_CCvalues[i]
		endif
	Endfor


	////////////// Plot //////////////
	duplicate/FREE/o List_CCvalues_All temp1
	Make/N=80/O Hist_All_ALL;DelayUpdate
	Histogram/B={-1,0.025,80} temp1,Hist_All_ALL

	duplicate/o List_Sign_CCvalues_All temp1
	Make/N=80/O Hist_Sign_All;DelayUpdate
	Histogram/B={-1,0.025,80} temp1,Hist_Sign_ALL

	Display/k=1 Hist_All_ALL
	appendtograph Hist_Sign_ALL
	ModifyGraph mode(Hist_All_ALL)=5,rgb(Hist_All_ALL)=(0,0,0)
	ModifyGraph mode=5,rgb=(0,0,0),hbFill(Hist_Sign_ALL)=2
	Legend/C/N=text1/J/A=RT "\\s(Hist_All_ALL) All correlation coefficient values\r\\s(Hist_Sign_All) Significant correlation coefficient values"

	SetDataFolder saveDFR

END

Function CorrSignSpearAuto(histograph,LWN,MWN,target)

	string histograph,LWN,MWN, target
	////// GUI ///////
	//prompt SPEARorPEAR_select, "SPEAR or PEAR?", popup "Spearman's Rank Correlation;Pearson's Correlation (assumes normally distributed data)"
	//prompt histograph, "Would you like me to plot a histogram of the correlation coefficient values?", popup "Yes please :);No, I'm good. Thanks though :)"
	//prompt LWN, "Select a 1D wave time-series", popup LocoWaveName
	//prompt MWN, "Select a 2D wave containing the calcium traces of your ROIs", popup MatrixWaveName
	//prompt Hz, "What is the acquisition frequency?"
	//doprompt "Correlation of time-series + significance test and mult. comp. corr.",SPEARorPEAR_select,histograph,LWN,MWN,Hz
	//	if(V_flag==1)
	//		abort
	//	endif
	string SPEARorPEAR_select = "Spearman's Rank Correlation"
	wave stimcycle
	variable sigpvalue = 0.20
	variable nshuffles = 1000
	variable pvalue = 0
	variable Hz=10

	print SPEARorPEAR_select
	variable SPEARorPEAR
	if (stringmatch(SPEARorPEAR_select,"Spearman's Rank Correlation"))
		SPEARorPEAR=0
	elseif (stringmatch(SPEARorPEAR_select,"Pearson's Correlation (assumes normally distributed data)"))
		SPEARorPEAR=1
	endif
	
	duplicate/FREE/o $MWN All
	duplicate/FREE/o $LWN Loco

	wave stim
	wave loco_timex

	// make a new folder "Correlation" if the folder does not aready exists
	string Correlation
	Correlation="root:" + target
	if (!datafolderexists(Correlation))
		newdatafolder $Correlation
	elseif (datafolderexists(Correlation))
		string 	message = "The folder '" + target + "' already exists"
		print message
	endif
	variable NumFrames=dimsize(Loco,0)
	make/FREE/o/n=(dimsize(All,1)) List_Corr_LOCO_CCvalues
	make/FREE/o/n=(dimsize(All,1)) List_Corr_LOCO_Pvalues

	Make/FREE/o/n=(NumFrames) RandomWave 
	Make/FREE/o/n=(NumFrames) RandomCorValue

	Make/FREE/o/n=0 ROI_list_Corr_NS
	Make/FREE/o/n=0 ROI_list_Corr_pos
	Make/FREE/o/n=0 ROI_list_Corr_neg

	variable LenRandom=dimsize(RandomWave,0)
	variable LenAll=dimsize(All,0)
	variable LenROIs=dimsize(All,1)

	// Generating a wave with a 1000 random starting points for the circular shifting
	RandomWave=0
	variable i,j,k,ii,jj,kk,iii,jjj,count1
	variable willekeur
	for (j=0;j<LenRandom;j+=1)
		i=0
		do
			willekeur=round(enoise(LenAll/2))+(LenAll/2)
			Findvalue/V=(willekeur) RandomWave
			If (V_Value==-1 && willekeur!=0)
				RandomWave[j]=willekeur
				i=1
			endif
		while (i!=1)
	endfor

	// Big Loop through all ROIs

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	// SPEARMAN                                                                                                      //
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	if (SPEARorPEAR==0)

		variable FirstValue
		variable HighValue
		//Correlation: Loop through all ROIs
		FOR (i=0;i<LenROIs;i+=1)
			count1=0
			variable CorrList
			duplicate/FREE/o/r=[][i] All temp1
			StatsRankCorrelationTest/ALPH=0.05/T=1/Q Loco, temp1
			wave W_StatsRankCorrelationTest
			CorrList=W_StatsRankCorrelationTest[4]
			List_Corr_LOCO_CCvalues[i]=CorrList
	
			//Correlation done with random shifts a number of times equal to RandomWave size
			for (j=0;j<LenRandom;j+=1)
				duplicate/FREE/o Loco tempStart
				Duplicate/FREE/o/r=[0,RandomWave[j]-1] Loco tempEnd
				Deletepoints/m=0 0,RandomWave[j]+0, tempStart 
				Concatenate/o/np=0 {tempStart,tempEnd}, tempRandom
				if (dimsize(tempRandom,0)!=dimsize(Loco,0)) // Wavelength mismatch
					duplicate/FREE/o Loco tempStart
					Duplicate/FREE/o/r=[0,RandomWave[j]-1] Loco tempEnd
					Deletepoints/m=0 0,RandomWave[j]+1, tempStart // Need to add +1 to make sure there is no wavelength mismatch
					Concatenate/o/np=0 {tempStart,tempEnd}, tempRandom
				endif
				StatsRankCorrelationTest/ALPH=0.05/T=1/Q tempRandom, temp1
				wave W_StatsRankCorrelationTest
				RandomCorValue[j]=W_StatsRankCorrelationTest[4] // 4 for Spearman				This wave now contains correlation values for nshuffle random shuffles
			endfor	
	
			//Calculate the distribution of corr values 
			Make/N=(LenRandom/4)/O Corr_Hist
			Histogram/cum/C/B=1 RandomCorValue, Corr_Hist 
			//Cumulative histogram
			Corr_Hist/=LenRandom										//Divide by number of samples to get CDF	
	
			if(CorrList<0)			
				pvalue = Corr_Hist(List_Corr_LOCO_CCvalues[i])		// Gives pvalue for negative corr values (analagous to the t-value): represents the area under the PDF curve to the LEFT of the corr value 
			else
				pvalue= 1 - Corr_Hist(List_Corr_LOCO_CCvalues[i])	// Take 1-CDF value at corr value (t-value) to acquire the area under the PDF (probability density function) to the RIGHT of the corr-value
			endif
	
			List_Corr_LOCO_Pvalues[i]=pvalue
		ENDFOR

		//stats_BHFDR(List_Corr_LOCO_Pvalues) // Correct for multiple comparison

		// Devide the ROIs in significantly positively or negatively correlated with locomotion, or not significant
		//Changed temporarily to make it more permissive. Original 0.05. Also change in the histogram.
		variable PvalueThr=0.05
		FOR (i=0;i<LenROIs;i+=1)
			if (List_Corr_LOCO_Pvalues[i]>PvalueThr)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_NS,0)),1, ROI_list_Corr_NS
				ROI_list_Corr_NS[dimsize(ROI_list_Corr_NS,0)-1]=i
			elseif (List_Corr_LOCO_Pvalues[i]<=PvalueThr && List_Corr_LOCO_CCvalues[i]>=0)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_pos,0)),1, ROI_list_Corr_pos
				ROI_list_Corr_pos[dimsize(ROI_list_Corr_pos,0)-1]=i
			elseif (List_Corr_LOCO_Pvalues[i]<=PvalueThr && List_Corr_LOCO_CCvalues[i]<0)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_neg,0)),1, ROI_list_Corr_neg
				ROI_list_Corr_neg[dimsize(ROI_list_Corr_neg,0)-1]=i
			endif
		ENDFOR

		string window_name=winname(0,2,1)
		killwindow $window_name
		//killwaves W_StatsRankCorrelationTest, tempRandom, Loco_Corrected

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
		// PEARSON                                                                                                         //
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
	elseif (SPEARorPEAR==1)

		FOR (i=0;i<LenROIs;i+=1)
			count1=0
			duplicate/FREE/o/r=[][i] All temp1
			make/FREE/o/n=1 TijdelijkeCorrList
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// Determine the highest correlation value by shifting the locomotion time-series minus to pls 2 seconds //
			duplicate/FREE/o Loco tempStart
			Duplicate/FREE/o/r=[(dimsize(Loco,0)-(2*Hz)),dimsize(Loco,0)-1] Loco tempEnd
			Deletepoints/m=0 (dimsize(Loco,0)-(2*Hz)),dimsize(Loco,0), tempStart
			Concatenate/o/np=0 {tempEnd,tempStart}, Loco_Corrected
			StatsLinearCorrelationTest/ALPH=0.05/Q Loco_Corrected, temp1
			wave W_StatsLinearCorrelationTest
			TijdelijkeCorrList[0]=W_StatsLinearCorrelationTest[1]
			for (ii=0;ii<4*Hz;ii+=1)
				FirstValue=Loco_Corrected[0]
				Deletepoints/m=0 0,1, Loco_Corrected
				insertpoints/m=0 dimsize(Loco_Corrected,0),1, Loco_Corrected
				Loco_Corrected[dimsize(Loco_Corrected,0)-1]=FirstValue
				StatsLinearCorrelationTest/ALPH=0.05/Q Loco_Corrected, temp1
				wave W_StatsLinearCorrelationTest
				insertpoints/m=0 dimsize(TijdelijkeCorrList,0),1, TijdelijkeCorrList
				TijdelijkeCorrList[dimsize(TijdelijkeCorrList,0)-1]=W_StatsLinearCorrelationTest[1]
			endfor
			wavestats/q TijdelijkeCorrList
			If (abs(V_max)>abs(V_min))
				HighValue=V_max
			else
				HighValue=V_min
			endif
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			StatsLinearCorrelationTest/ALPH=0.05/Q Loco, temp1
			wave W_StatsLinearCorrelationTest
			List_Corr_LOCO_CCvalues[i]=W_StatsLinearCorrelationTest[1] // 1 for Pearson
			for (j=0;j<nshuffles;j+=1)
				duplicate/FREE/o Loco tempStart
				Duplicate/FREE/o/r=[0,RandomWave[j]-1] Loco tempEnd
				Deletepoints/m=0 0,RandomWave[j], tempStart
				Concatenate/o/np=0 {tempStart,tempEnd}, tempRandom
				if (dimsize(tempRandom,0)!=dimsize(Loco,0)) // Wavelength mismatch
					duplicate/FREE/o Loco tempStart
					Duplicate/FREE/o/r=[0,RandomWave[j]-1] Loco tempEnd
					Deletepoints/m=0 0,RandomWave[j]+1, tempStart // Need to add +1 to make sure there is no wavelength mismatch
					Concatenate/o/np=0 {tempStart,tempEnd}, tempRandom
				endif
				StatsLinearCorrelationTest/ALPH=0.05/T=1/Q tempRandom, temp1
				wave W_StatsLinearCorrelationTest
				RandomCorValue[j]=W_StatsLinearCorrelationTest[1] // 1 for Pearson
			endfor
	
			//Calculate the distribution of corr values 
			Make/N=(LenRandom/4)/O Corr_Hist
			Histogram/CUM/C/B=1 RandomCorValue, Corr_Hist     //Cumulative histogram
			Corr_Hist/=LenRandom
			//Divide by number of samples to get CDF	
			if(CorrList<0)			
				pvalue = Corr_Hist(List_Corr_LOCO_CCvalues[i])		// Gives pvalue for negative corr values (analagous to the t-value): represents the area under the PDF curve to the LEFT of the corr value 
			else
				pvalue= 1 - Corr_Hist(List_Corr_LOCO_CCvalues[i])	// Take 1-CDF value at corr value (t-value) to acquire the area under the PDF (probability density function) to the RIGHT of the corr-value
			endif
	
			List_Corr_LOCO_Pvalues[i]=pvalue
		ENDFOR

		//stats_BHFDR(List_Corr_LOCO_Pvalues) // Correct for multiple comparison

		// Devide the ROIs in significantly positively or negatively correlated with locomotion, or not significant
		FOR (i=0;i<LenROIs;i+=1)
			if (List_Corr_LOCO_Pvalues[i]>0.05)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_NS,0)),1, ROI_list_Corr_NS
				ROI_list_Corr_NS[dimsize(ROI_list_Corr_NS,0)-1]=i
			elseif (List_Corr_LOCO_Pvalues[i]<=0.05 && List_Corr_LOCO_CCvalues[i]>=0)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_pos,0)),1, ROI_list_Corr_pos
				ROI_list_Corr_pos[dimsize(ROI_list_Corr_pos,0)-1]=i
			elseif (List_Corr_LOCO_Pvalues[i]<=0.05 && List_Corr_LOCO_CCvalues[i]<0)
				InsertPoints/M=0 (dimsize(ROI_list_Corr_neg,0)),1, ROI_list_Corr_neg
				ROI_list_Corr_neg[dimsize(ROI_list_Corr_neg,0)-1]=i
			endif
		ENDFOR

		string window_name2=winname(0,2,1)
		killwindow $window_name2
		killwaves W_StatsLinearCorrelationTest, tempRandom, Loco_Corrected

	endif
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Put the new waves in the "Correlation" Folder

	string	stringN = Correlation + ":ROI_list_Corr_neg"
	duplicate/o ROI_list_Corr_neg $stringN
	stringN= Correlation +  ":ROI_list_Corr_pos"
	duplicate/o ROI_list_Corr_pos $stringN
	stringN=	Correlation + ":ROI_list_Corr_NS"
	duplicate/o ROI_list_Corr_NS $stringN

	stringN = Correlation + ":List_Corr_LOCO_CCvalues"
	duplicate/o List_Corr_LOCO_CCvalues $stringN
	stringN = Correlation + ":List_Corr_LOCO_Pvalues"
	duplicate/o List_Corr_LOCO_Pvalues $stringN
	stringN=Correlation + ":Locomotion"
	duplicate/o Loco $stringN
	stringN = Correlation + ":All_CorrectedQA_DivMoff"
	duplicate/o All $stringN

	if (stringmatch(histograph,"Yes"))
		TwoClusterHisttarget(target)
	endif

END


//Function TwoClusterHist()
//
//wave Root:Correlation:List_Corr_LOCO_CCvalues
//wave Root:Correlation:List_Corr_LOCO_Pvalues
//wave Root:Correlation:ROI_list_Corr_neg
//wave Root:Correlation:ROI_list_Corr_pos
//
//Make/o/n=(dimsize(Root:Correlation:List_Corr_LOCO_CCvalues,0))  Root:Correlation:List_CCvalues_All
//Make/o/n=0  Root:Correlation:List_Sign_CCvalues_All
//variable LenListCC=dimsize(Root:Correlation:List_Corr_LOCO_CCvalues,0)
//
//
//variable i,ii,j,jj,k,kk
//
//For (i=0;i<LenListCC;i+=1)
//	 Root:Correlation:List_CCvalues_All[i]=Root:Correlation:List_Corr_LOCO_CCvalues[i]
//	if (Root:Correlation:List_Corr_LOCO_Pvalues[i]<0.05)
//		insertpoints/m=0 dimsize(Root:Correlation:List_Sign_CCvalues_All,0),1, Root:Correlation:List_Sign_CCvalues_All
//		Root:Correlation:List_Sign_CCvalues_All[dimsize(Root:Correlation:List_Sign_CCvalues_All,0)-1]=Root:Correlation:List_Corr_LOCO_CCvalues[i]
//	endif
//Endfor
//
//
//////////////// Plot //////////////
//duplicate/o Root:Correlation:List_CCvalues_All temp1
//Make/N=80/O Root:Correlation:Hist_All_ALL;DelayUpdate
//Histogram/B={-1,0.025,80} temp1,Root:Correlation:Hist_All_ALL
//
//duplicate/o Root:Correlation:List_Sign_CCvalues_All temp1
//Make/N=80/O Root:Correlation:Hist_Sign_All;DelayUpdate
//Histogram/B={-1,0.025,80} temp1,Root:Correlation:Hist_Sign_ALL
//
//Display/k=1 Root:Correlation:Hist_All_ALL
//appendtograph Root:Correlation:Hist_Sign_ALL
//ModifyGraph mode(Root:Correlation:Hist_All_ALL)=5,rgb(Root:Correlation:Hist_All_ALL)=(0,0,0)
//ModifyGraph mode=5,rgb=(0,0,0),hbFill(Root:Correlation:Hist_Sign_ALL)=2
//TextBox/C/N=text0/A=RT "Cluster #1"
//
//
//END

function TR_CorrectAndCorrelateFiles()

	string objList="CroppedCh1_DFF0;zestimate;zestimateflipped;locotrace;stimwave;CroppedRef_Ch1_reg_SD_QA;zestimate;"
	string pathstr="C:\Data\xCorrect\VIP_Boutons\ToCorrelate"
	newpath/O path, pathstr

	string filelist= indexedfile(path,-1,".pxp")	
	filelist=SortList(filelist, ";", 16)				  	  	//Create a list of all the .pxp files in the folder. -1 parameter addresses all files.
		
	Variable Length=itemsinlist(filelist),i
	print pathstr
		
	for(i=0;i<Length;i+=1)	// cycle through files in the folder
		
		string expName= StringFromList(i, filelist)
		string message = "Loading from Experiment " + expName + " ( " + num2str(i) + "/" +num2str(length) + " in this folder"
		print message
		
		string	expPath = pathstr + "\\" + expName
		LoadData/J=objList   /Q expPath
		
		
		string LON="LocoTrace"
		string STN="stimWave"
		string locoN=wavelist("Rec*_Ch3","","")
		string stimN=wavelist("Rec*_Ch4","","")
		string zest=wavelist("zestimate","","")
		string MWN="root:correctedQA_DivMoff"

		correctqa()

		CorrSignSpearAuto("No",LON,MWN,"Correlation_CorrectedCurated")
		//TextBox/C/N=text0/A=MC "Corrected Curated"
	
		MWN="root:pre_CorrectedQA"
		CorrSignSpearAuto("No",LON,MWN,"Correlation_RawCurated")
		//	TextBox/C/N=text0/A=MC "Raw Curated"
	
		MWN="root:CroppedCh1_DFF0"
		CorrSignSpearAuto("No",LON,MWN,"Correlation_RawAll")
		//	TextBox/C/N=text0/A=MC "Raw All"
	
		MWN="root:All_CorrectedQA_DivMoff"
		CorrSignSpearAuto("No",LON,MWN,"Correlation_CorrectedAll")
		//	TextBox/C/N=text0/A=MC "Corrected All"
		duplicate/o $zest zestW

		zestW*=-1
		duplicate/o zestW $zest
	
	
		correctqa()
		zestW*=-1
		duplicate/o zestW $zest
		zestW*=-1

		MWN="root:correctedQA_DivMoff"
		CorrSignSpearAuto("No",LON,MWN,"Correlation_CorrectedCurated_inv")
		//TextBox/C/N=text0/A=MC "Corrected Curated_inv"
	
		MWN="root:pre_CorrectedQA"
		CorrSignSpearAuto("No",LON,MWN,"Correlation_RawCurated_inv")
		//	TextBox/C/N=text0/A=MC "Raw Curated_inv"
	
		MWN="root:All_CorrectedQA_DivMoff"
		CorrSignSpearAuto("No",LON,MWN,"Correlation_CorrectedAll_inv")
		//	TextBox/C/N=text0/A=MC "Corrected All_inv"	
	
	
			
		string saveName = expName
		message = "Saving experiment"
		saveExperiment/P=path as saveName
		killallwaves()
	endfor
end



function TR_CorrectAndCorrelate()

	string LON="locotrace"
	string zest=wavelist("zestimate","","")
	string MWN="root:correctedQA_DivMoff"

	correctqa()

	CorrSignSpearAuto("Yes",LON,MWN,"Zest_Correlation_CorrectedCurated")
	TextBox/C/N=text0/A=MC "Corrected Curated"
	
	MWN="root:pre_CorrectedQA"
	CorrSignSpearAuto("Yes",LON,MWN,"Zest_Correlation_RawCurated")
	TextBox/C/N=text0/A=MC "Raw Curated"
	
	MWN="root:CroppedCh1_SD_QA"
	CorrSignSpearAuto("Yes",LON,MWN,"Zest_Correlation_RawAll")
	TextBox/C/N=text0/A=MC "Raw All"
	
	MWN="root:All_CorrectedQA_DivMoff"
	CorrSignSpearAuto("Yes",LON,MWN,"Zest_Correlation_CorrectedAll")
	TextBox/C/N=text0/A=MC "Corrected All"

	
end


Function TwoClusterHistPooled(CValues,PValues)

wave CValues, PValues
duplicate/o CValues List_Corr_LOCO_CCvalues
duplicate/o PValues List_Corr_LOCO_Pvalues
	
	variable LenListCC=dimsize(List_Corr_LOCO_CCvalues,0)
	
	Make/FREE/o/n=(LenListCC)  List_CCvalues_All
	Make/o/n=0  List_Sign_CCvalues_All

	variable i
	
	variable PvalueThr=0.05
	
	For (i=0;i<LenListCC;i+=1)
		List_CCvalues_All[i]=List_Corr_LOCO_CCvalues[i]
		if (List_Corr_LOCO_Pvalues[i]<PvalueThr)
			insertpoints/m=0 dimsize(List_Sign_CCvalues_All,0),1, List_Sign_CCvalues_All
			List_Sign_CCvalues_All[dimsize(List_Sign_CCvalues_All,0)-1]=List_Corr_LOCO_CCvalues[i]
		endif
	Endfor


	////////////// Plot //////////////
	duplicate/FREE/o List_CCvalues_All temp1
	Make/N=80/O Hist_All_ALL;DelayUpdate
	Histogram/B={-1,0.025,80} temp1,Hist_All_ALL

	duplicate/o List_Sign_CCvalues_All temp1
	Make/N=80/O Hist_Sign_All;DelayUpdate
	Histogram/B={-1,0.025,80} temp1,Hist_Sign_ALL

	Display/k=1 Hist_All_ALL
	appendtograph Hist_Sign_ALL
	ModifyGraph mode(Hist_All_ALL)=5,rgb(Hist_All_ALL)=(0,0,0)
	ModifyGraph mode=5,rgb=(0,0,0),hbFill(Hist_Sign_ALL)=2
	//Legend/C/N=text1/J/A=RT "\\s(Hist_All_ALL) All correlation coefficient values\r\\s(Hist_Sign_All) Significant correlation coefficient values"
	
	SetAxis bottom -0.6,0.6
SetAxis left *,80
Label left "No. ROIs"
ModifyGraph fSize=14,axThick=1.5;DelayUpdate
Label left "\\Z16\\f01No. ROIs";DelayUpdate
Label bottom "\\Z16\\f01Correlation Coefficient"
ModifyGraph axOffset(left)=2,axOffset(bottom)=-0.5,standoff=0

ModifyGraph axOffset(bottom)=-0,lblMargin(left)=20
	
	
	// compute how many positive, negative and NS ROIs there are and print them for the user
	variable numSigPos=0,numSigNeg=0
	for(i=0;i<(dimsize(List_Sign_CCValues_All,0));i+=1)
		if(List_Sign_CCValues_All[i]>0)
			numSigPos+=1
		else
			numSigNeg+=1
		endif
	endfor
	
	print "There are " + num2str(numSigNeg) + " significantly negatively correlated ROIs"
	
	print "There are " + num2str(numSigPos) + " significantly positively correlated ROIs"

END
