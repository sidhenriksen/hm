# Response of V1 cells to half-matched random dot stereograms

## The browser
This code includes an interactive data browser that lets you explore a curated
version of the dataset. However, because Matlab has made changes to how figures
and callback functions are handled in recent versions, the data browser is very
dodgy for Matlab versions earlier than 2014b. If you have an older version of Matlab,
please see the "The data" section for details about accessing the raw data itself.

In order to use the browser, clone/download the repository and run
```
PlotCurated
```
which will open an interactive figure with two subplots. The left subplot shows data for 5% density, while
the right subplot shows for 24% density. Each point in the plot shows a single cell. By default, the axes 
will show the r-value (correlation coefficient) between the anticorrelated and correlated tuning curves on 
the y-axis, and between the half-matched and correlated tuning curves on the x-axis. You can change the
axes by selecting "Plot -> Main: x axis", and choosing an appropriate metric (for example half-matched slope).
Most metrics should be self-explanatory upon reading Henriksen, Read, & Cumming (2016), but should you have
any questions about the data, please get in touch [sid(dot)henriksen(at)gmail(dot)com].

In the plots, each unique colour shows a different recording session. The size of the points shows the magnitude of
the Disparity Discrimination Index. Square points correspond to cells from monkey Jbe, while circular points correspond 
to cells from monkey Lem. Clicking the dots will bring up tuning curves showing that cells' responses to correlated (red line),
anticorrelated (black line) and half-matched stereograms (blue line). The disparity tuning curve plot also allows you to
toggle different error bars (SEM, SD, and 95% bootstrap CIs). 

## The data
The data browser requires Matlab version 2014b or later. If you don't have access to this (or you don't have Matlab at all),
you can access the curated data itself and explore it manually. The data is stored in the file `CuratedCells.mat`. If you run
```
load('CuratedCells.mat')
```
while in the repo directory, a variable called `Base` will be loaded into your workspace. `Base` is a struct with fields
`density`, `Cells`, `exptlist`, `penetrationlist`. Base is a 1x2 array, where the elements correspond to the two densities
used (5% and 24%).

`density` tells the density of the random dot stereograms used (5% or 24%).
`Cells` is the struct containing all the curated data.
`exptlist` lists all the recording sessions
`penetrationlist` gives the x-y coordinates for each recording session.

### The Cells struct
This contains all the curated data for these experiments. To access this data:
```
density5_cell1 = Base(1).Cells(3);
```
This gives you the data for the 1st cell, where the dot density was 5%. Base(1)
corresponds to 5% density recordings, and Base(2) to 24% density recordings.
To plot a tuning curve we can do:
```
current_cell = Base(1).Cells(1);
dx = current_cell.Dx;
correlated = current_cell.correlatedResponse;
anticorrelated = current_cell.anticorrelatedResponse;
halfmatched = current_cell.halfmatchedResponse;
figure();
plot(dx,correlated,'r -',dx,halfmatched,'- b',dx,anticorrelated,'k -','linewidth',3);
xlabel('Disparity (deg)');
ylabel('Spike count');
```

The following gives a complete documentation of the fields in Base(k).Cells:
`cellnumber` - the cell number given in the recording session.  
`filename` - where the data is located locally (not useful for external use).  
`regHm` - output of a type 2 regression between correlated and half-matched tuning curves.
This is a vector of size three, where the entries are [r,m,b]. r is the correlation coefficient,
m is the half-matched slope, b is the offset (probably useless).  
`regAc` - same as regHm just for anticorrelated slope.  
`regHmRegular` - same as regHm except using OLS regression.  
`regAcRegular` same as regAc except using OLS regression.  
`Dx` - disparities used in the experiment.

`correlatedResponse` - tuning curve for correlated stimuli (averaged across trial)  
`halfmatchedResponse` - tuning curve for half-matched stimuli  
`anticorrelatedResponse` - tuning curve for anticorrelated stimuli  

`correlatedSEM` - SEMs for the correlated tuning curve  
`halfmatchedSEM` - SEMs for the half-matched tuning curve  
`anticorrelatedSEM` - SEMs for the anticorrelated tuning curve  

`RMS` - root mean square for tuning curves for correlated, half-matched, and anticorrelated stimuli  
`DDI` - DDI computed on correlated tuning curves  
`HMauc` - area under the ROC curve (AUROC) for half-matched stimuli  
`HMdprime` - HMauc converted to a d' value  
`Cauc` - AUROC for correlated stimuli  
`Cdprime` - Cauc converted to a d' value  
`HMaucZ` - AUROC for half-matched stimuli (spike counts Z-normalised by block number to correct for slow drifts)  
`CaucZ` - same as `HMaucZ` for correlated stimuli  

`density` - density of the stimulus  
`dw` - dot width used (in degrees)  
`ciLowHm` and `ciHighHm` - 95% bootstrap confidence intervals (CIs) for half-matched tuning curves  
`ciLowC` and `ciHighC` - 95% bootstrap CIs for correlated tuning curves  
`ciLowAc` and `ciHighAc` 95% bootstrap CIs for anticorrelated tuning curves  
`CHm_r_CI` - 95% bootstrap CIs for r between correlated and half-matched  
`CHm_slope_CI` - 95% bootstrap CIs for half-matched slope  
`CAc_r_CI` - 95% bootstrap CIs for r between correlated and anticorrelated  
`CAc_slope_CI` - 95% bootstrap CIs for anticorrelated slope  
`DDIhm` - half-matched DDI  