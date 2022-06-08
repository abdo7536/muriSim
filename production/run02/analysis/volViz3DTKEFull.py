################# Script to load LegacyVTK/PVTU visualization files into the pipeline and conduct Volume rendering on the user defined field(s)
################# Input(s): Sorted by python-paraview custom library functions. NOTE: All inputs listed in this section are mandatory
########## General Input(s):
filnam = "/p/work1/abdo7536/SHIT01/subvol1*TKE.vtk"
### extSS = Switch to enable extraction of subset of data (1 - Yes; 0 - No); NOTE: If set to 1 then prescribe inputs for 'extSubset' funciton below
extSS = 0 # Extract subset from the subvol_xxxx.vtk domain
### transForm = Switch to enable data transformation to the visualization data (1 - Yes; 0 - No); NOTE: If set to 1 then prescribe inputs for 'pyCalc' function below
transForm = 0 # Transform visualization data 
### datRng = The data range for Visualization data field - ordered [min, max]
datRng = [0.0, 0.1]
### invColorTF = Switch that inverts the color transfer function (1 -- invert; 0 - dont invert)
invColorTF = 0 # Color TF invert
### opcFn = The opacity transfer function option: 
### 1 - Hyperbolic Tangent; 2 - Linearly decreasing; 3 - Linearly increasing; 4. Single step (Increasing); 5 - Four steps (Increasing); 
### 6 - Wedge; 7 - Gaussian; 8 - Four steps (Increasing); 0 - Custom opacity TF weights 
### NOTE: Both functions decrease in opacity Wt. with increasing data values
opcFn = 3
### opcTbl = Input used to prescribe custom opacity table NOTE: Use this only if opcFn = 0 was set above
#opcTbl = [-0.2, 1.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0]
### opcPts = The number of points on the Opacity Transfer Function Between the data ranges (for ranges see the PDF of data)
opcPts = 100
### scalOpc = Scaling Factor applied to the Transfer Function: 0 < scalOpc < 1 (NOTE: This sets the maximum value of Opacity applied to the colormap)
scalOpc = 0.75
### modColorTF = Switch to enable modifications to the default color TF; NOTE: This is an advanced feature, recommended that a general user sets this input to 0
# NOTE: If set to 1 then ColorTF is a mandatory input
modColorTF = 0 # Modify the color transfer function
### ColorTF = input containing modified Color Transfer function; NOTE: Arranged as [datapoint1, R1, G1, B1, datapoint2, R2, G2, B2, ....]
### EXAMPLE INPUT: ColorTF [0.0, 0.231373, 0.298039, 0.752941, 4.0, 0.865003, 0.865003, 0.865003, 6.479750633239746, 0.9490196078431372, 0.5647058823529412, 0.4470588235294118, 8.0, 0.705882, 0.0156863, 0.14902]
### renVieBG = The background colour for the render view: Set as RGB colour weights [0.0, 0.0, 0.0] - for Black and [1.0, 1.0, 1.0] - for white
renVieBG = [0.0, 0.0, 0.0]
### hor_res = Horizontail Pixel Count for the desired render resolution # for 4K res: 3840 pixels; for 8K res: 7,680 pixels [type - int]
hor_res = 7680
### vert_res = Vertical Pixel Count for the desired render resolution # for 4K res: 2160 pixels; for 8K res: 4,320 pixels [type - int]
vert_res = 4320
# txtAnntn = Switch to enable text annotation (1 - on; 0 - off); NOTE: Can be used to print time annotation or simulation details
txtAnntn = 0    # EXAMPLE INPUT

########## function name: 'lodatsin'
### fname = Name of the file to be loaded into Paraview pipeline [type - string]
# fname ----- declared in the program below 
### colrby = The data field name as stored in the vtk file [type - string]
colrby = 'KE-diss-rate'         # [common options are: 'KE-diss-rate','Lambda2','T-diss-rate','p\'','u\'','v\'','w\'','vort_x','vort_y','vort_z','vort-mag','stab-index']
### rv = It is the 'renderView' of the paraview layout in which the visualization is being carried out [type - paraview object]
# rv ----- Obtained from the step where paraview layout and render view are set
### colschem = Set the user desired colourmap -- currently only works for paraview presets [type - string]
### Some examples include: Linear: 'Inferno (matplotlib)', 'Viridis (matplotlib)', 'Magma (matplotlib)', 'Plasma (matplotlib)'
### Default matlab schemes: 'Jet', 'Jet', 'Cool to Warm', 'Cool to Warm (Extended)', 'Rainbow Desaturated', 'Rainbow Uniform'
colschem = 'GnYlRd'
### form = File extension [type - string]
form = 'vtk' # [EXAMPLE INPUT]. NOTE: Only two possible inputs - "vtk" or "pvtu"
### lodFlds = Field names contained in the PVTU file: Choose from the list of:
# 'KE-diss-rate','Lambda2','T-diss-rate','p\'','T\'','u\'','v\'','w\'','vort_x','vort_y','vort_z','vort-mag','stab-index'
lodFlds = ['KE-diss-rate']

########## function name: 'extSubset'
### domXYZ = The Simulation Box size (NOTE: See inputs for Xl, Yl, Zl, and Zo in sam.inp) [NOTE: list of integers]
domXYZ = [120, 42, 24]
### x = The x1-x2 Extent for the plane/volume (NOTE: The domain spans -X/2 --> X/2; x1 and x2 are +ve distances from -X/2) [type - list of integers]
x = [35, 70]
### y = The y1-y2 Extent for the plane/volume (NOTE: The domain spans -Y/2 --> Y/2; y1 and y2 are +ve distances from -Y/2) [type - list of integers]
y = [15, 42]
### z = The z1-z2 Extent for the plane/volume (NOTE: The domain spans -Z/2 --> Z/2; z1 and z2 are +ve distances from Zo) [type - list of integers]
z = [6, 18]
### N = The number of mesh points in X, Y, and Z for the simulation Domain (NOTE: Inputs Nx, Ny, and Nz - defined in sam.inp file) [type: list of integer$
N = [7200, 2100, 1200]
### colrby = The field by which the plane-data is coloured [type - string]
# colrby ----- see input for function 'lodatsin'
### objnam = The pipeline browser object on which the subset operation is being carried out on [type - paraview object]
# objnam ----- see input for function 'lodatsin'
### rv1 = Render view object [type - paraview object]
# rv ----- Obtained from the step where paraview layout and render view are set
### subSamp = The subsampling rate applied to the data [type - list of integers]. NOTE: [1, 1, 1] - for no subsampling (Default: [1, 1, 1])
subSamp = [1, 1, 1]

########## function name: 'pyCalc'
### objnam = The pipeline browser object on which the subset operation is being carried out on [type - paraview object]
# objnam ----- see input for function 'lodatsin'
### paramNm = Name of the calculated parameter/variable [type - string] (use single quotes to define the expression)
paramNm = colrby
### expr = The calculation expression; NOTE: All the variables used here should be available in the 'objnam' data object (use single quotes to define the exp$
expr = '(Lambda2*-1)**0.5'
### rv1 = Render view object [type - paraview object]
# rv ----- Obtained from the step where paraview layout and render view are set

########## function name: 'modclrbar'
### Input(s): Mandatory
### RV1 = Render View object [type - paraview object]
# RV1 ----- Passed to the function in the main program
### colrby = The data field name with which the data ought to be coloured by [type - string]
# colrby ----- Already defined above
### Input(s): Optional
### FontSize = Text and label font size [type - int] (DEFAULT: 4)
FontSize = 4
### ornt = The orientation of the colorbar [type - string] (DEFAULT: 'Horizontal')
ornt = 'Vertical'             # 'Horizontal' or 'Vertical'
### sclPts = Number of ticks to display on colorbar
sclPts = 5
### cbarpos_H = Horizontal position of the colorbar [type - int/float] (DEFAULT: 0.5)
cbarpos_H = 0.025
### cbarpos_V = Vertical position of the colorbar [type - int/float] (DEFAULT: 0.01)
cbarpos_V = 0.4
### cbar_thic = Colorbar thickness [type - int/float] (DEFAULT: 3)
cbar_thic = 3
### cbar_len = Colorbar length [type - int/float] (DEFAULT: 0.25)
cbar_len = 0.25
### schem = The Paraview color scheme preset [type - string] (DEFAULT: 'Jet')
# schem ----- Already defined above (see input 'colschem')

########## function name: 'set_cam' 
### RV1 = Render View object [type - paraview object] 
# RV1 ----- Passed to the function in the main program
### rotCam = Input to control Camera Rotations (if desired); 
# 0 - No rotations
# 1 - Rotate camera about Roll/Azimuth/Elevation Axes for user defined Starting Camera position 
# NOTE: The Camera starting position in case of rotCam = 2 or 3 is XZ plane with +Z pointing up and +Y pointing out of the display
rotCam = 0

### If rotCam = 0; plnam, ptup, pos_ax are Mandatory but set other inputs as described below 
### plnam = 1 or 2 or 3 -- 1=xy planes are being visualized; 2=yz planes are being visualized; 3=xz planes are being visualized [type - int]
### ptup = The axis which is pointing up -- 1,-1 or 2,-2 or 3,-3; 1,-1 = X or -X axis; 2, -2 = Y or -Y axis; 3, -3 = Z or -Z axis [type - int] 
### pos_ax = The axis which is pointing into the display -- 1,-1 or 2,-2 or 3,-3; 1,-1 = X or -X axis; 2, -2 = Y or -Y axis; 3, -3 = Z or -Z axis [type - int] 
### rotax = Axis/Axes about which the camera rotates [NOTE: Can choose multiple]; A 1D array of 6 mandatory inputs (either 0 or 1 to turn on and off) arranged as [+X, +Y, +Z, -X, -Y, -Z] [type - listof int]
# rotax ----- set to [0,0,0,0,0,0]
### angdiv = The number of divisions between 0 - 360 degrees for each axes rotation; A 1D array of 6 integer values (preferably divisors/multiples of 360) [type - list of ints] 
# angdiv ----- set to 0
### savnam = Path to the save output .png files [type - listof int]
# savnam ----- set to "default"

### If rotCam = 1; stCampt, rotax, angdiv, savnam are important but set other inputs as described below 
### stCampt = The starting position of the camera (prescribed as Azimuth from +Y axis: -ve value for CCW rotations; Elevation from +Y axis: +ve for elevate up) [type - list of floats]
# NOTE: Arranged as [azm, elev, zoom] NOTE: set zoom value as 1 if no zoom is needed; >1 for zoomIn and <1 for zoomOut
stCampt = [-45, 90.0, 1.0]
### rotax = Input that prescribes the camera movement for roll, azimuth, elevation axes respectively, and camera zoom [NOTE: Can choose multiple inputs simultaneously]; A 1D array of 4 mandatory inputs (either 0 or 1 to turn on and off) arranged as [roll, azm, elev, zoom] [type - list of int]
rotax = [0, 1, 0, 0]
### angdiv = The number of divisions between 0 - 360 degrees for each axes rotation [type - integer] 
angdiv = [6,0,6]
### savnam = Path to the save output .pvtu files [type - listof int]
savrotnam = '/p/work1/abdo7536/SHIT01/TKEFull.'
### plnam = 1 or 2 or 3 -- 1=xy planes are being visualized; 2=yz planes are being visualized; 3=xz planes are being visualized [type - int]
# plnam ----- set to 0 
plnam = 1
### ptup = The axis which is pointing up -- 1,-1 or 2,-2 or 3,-3; 1,-1 = X or -X axis; 2, -2 = Y or -Y axis; 3, -3 = Z or -Z axis [type - int] 
# ptup ----- set to 0
ptup = 2
### pos_ax = The axis which is pointing into the display -- 1,-1 or 2,-2 or 3,-3; 1,-1 = X or -X axis; 2, -2 = Y or -Y axis; 3, -3 = Z or -Z axis [type - int] 
# pos_ax ----- set to 0
pos_ax = -3
### hor_res = Horizontail Pixel Count for the desired render resolution # for 4K res: 3840 pixels; for 8K res: 7,680 pixels [type - int] 
# hor_res ----- Already defined above
### vert_res = Vertical Pixel Count for the desired render resolution # for 4K res: 2160 pixels; for 8K res: 4,320 pixels [type - int] 
# vert_res ----- Already defined above
### stTranpBG =	Save images with transparent render view background: 0 - NO; 1 - YES [type - int]
stTranpBG = 0   # (if set to 0 then define in general inputs what the background should be use input 'renVieBG' and specify input as RGB weights)

########## function name: 'addtxt'
### Input(s): Mandatory
### rv1 = Render View object (a paraview.simple object which renderes the layout onscreen display) [type - paraview object]
# rv1 ----- Passed to the function after setting the paraview layout and render view 
### txt_str = The text to be displayed (By default this will be displayed at the upper right corner) [type - string]
txt_str = r't/$\tau_{b} =$ '
### txt_ftsz = The text source font size [type - int]
txt_ftsz = 120
### txt_loc = The location of annotation text 'LowerRightCorner', 'UpperRightCorner', 'UpperCenter', 'UpperLeftCorner', 'LowerLeftCorner', 'LowerCenter'
txt_loc = 'UpperLeftCorner'

################# import all the python-paraview essential functions (custom library)
from pvessentials import *

################# PROGRAM START
################# List the filenames to process
fls = []
fls = sorted(glob.glob(filnam))
################# Grab the t/Tb values from file
################# read data from 'x_spec.xxxx' file
if txtAnntn == 1:
    with open('tStTb01.txt') as etst:
        temparr = etst.readlines()
    tmpdat = np.zeros(len(temparr)-1)
    ### dump data from file into data arrays
    for ind in range(1,len(temparr)):
        tTb = temparr[ind]
        tmpdat[ind-1] = float(tTb[20:30])

################# Run a loop over the list of all files
for lp in range(0,len(fls)):
    ################# Step 0: Setup paraview camera, render view, and layout
    ### Processing time calculator - sample 0
    timer0 = time.perf_counter()
    ### Disable render on startup [Saves lot of time to not render every change]
    paraview.simple._DisableFirstRenderCameraReset()
    ### get active view
    RV1 = GetActiveViewOrCreate('RenderView')
    ### Set Background Color
    RV1.Background = renVieBG
    ### get the basic layout
    layout1 = GetLayout()    
    ### set active view
    SetActiveView(RV1)
    ### Processing time calculator - sample 1
    timer1 = time.perf_counter()
    print(f'Set Paraview Layout in {round(timer1-timer0, 2)} s')

    ################# Step 1: Load the plane(s) data into the pipeline
    [subSetobj,subSetobjDisp] = lodatsin(fls[lp],colrby,RV1,colschem,form,lodFlds)
    ### update the view to ensure updated data information
    RV1.Update()
    ### Processing time calculator - sample 2
    timer2 = time.perf_counter()
    print(f'Loaded file into the pipeline in {round(timer2-timer1, 2)} s')

    ################# Step (Optional): Extract subset volume
    if extSS == 1:
        [vol,volDisp] = extSubset(domXYZ,x,y,z,N,subSetobj,RV1,colrby,subSamp)
        RV1.Update()
    else:
        vol = subSetobj
        volDisp = subSetobjDisp
        RV1.Update()

    ################# Step 2 (optional): Transform the Lambda2 field using the Calculation filter
    if (transForm == 1):
        [calcObj,calcObjDisp] = pvCalc(vol,paramNm,expr,RV1)
        ### update the view to ensure updated data information
        RV1.Update()
        ### Processing time calculator - sample 3
        timer3 = time.perf_counter()
        print(f'Transformed L2 dataset in {round(timer3-timer2, 2)} s')
    else:
        ### Processing time calculator - sample 3
        timer3 = time.perf_counter()

    ################# Step 3: Change Representation to Lambda2 and Representation Type to Volume
    ### Change representation type
    if transForm == 1:
        calcObjDisp.SetRepresentationType('Volume')
    else:
        volDisp.SetRepresentationType('Volume')
    RV1.EnableRayTracing = 1
    ### Show AxesGrid
    RV1.AxesGrid.Visibility = 1
    RV1.AxesGrid.XTitle = ''
    RV1.AxesGrid.YTitle = ''
    RV1.AxesGrid.ZTitle = ''
    RV1.AxesGrid.ShowTicks = 0
    RV1.AxesGrid.CullBackface = 1
    RV1.AxesGrid.CullFrontface = 1
    ### get color legend/bar for the field in RV1
    if transForm == 1:
        [TF, TFbar] = modclrbar(RV1, paramNm, FontSize, ornt, cbarpos_H, cbarpos_V, cbar_thic, cbar_len, colschem)
        ### get Colour Transfer Function
        colTF = GetColorTransferFunction(paramNm)
    else:
        [TF, TFbar] = modclrbar(RV1, colrby, FontSize, ornt, cbarpos_H, cbarpos_V, cbar_thic, cbar_len, colschem)
        ### get Colour Transfer Function
        colTF = GetColorTransferFunction(colrby)
    if invColorTF == 1:
        ### invert the Colour transfer function
        colTF.InvertTransferFunction()
    if transForm == 1:
        ### Get opacity transfer function/opacity map 
        TFopc = GetOpacityTransferFunction(paramNm)
        ### show color bar/color legend
        calcObjDisp.SetScalarBarVisibility(RV1, True)
    else:
        ### Get opacity transfer function/opacity map 
        TFopc = GetOpacityTransferFunction(colrby)
        ### show color bar/color legend
        #volDisp.SetScalarBarVisibility(RV1, True)

    ################# rescale color map to fit the current data range
    TF.RescaleTransferFunction(datRng[0], datRng[1])          # modifies the Colorbar TF
    TFopc.RescaleTransferFunction(datRng[0], datRng[1])       # modifies the opacity TF

    ################# Apply custom colorbar labels
    TFbar.UseCustomLabels = 1
    TFbar.CustomLabels = np.linspace(datRng[0],datRng[1],sclPts)
    ### update the view to ensure updated data information
    RV1.Update()
    ### Processing time calculator - sample 4
    timer4 = time.perf_counter()
    print(f'Applied Volume Representation in {round(timer4-timer3, 2)} s')

    ################# Step 4: Modify the Opacity Transfer Function
    ### Check if a custom opacity table was requested, if not apply the preset
    if opcFn!=0:
        ### Obtain the points on the opacity Transfer Function
        opcTFpts = np.linspace(datRng[0],datRng[1], num=opcPts)
        ### Calculate the opacity transfer function values/weights (NOTE: Different Opacity models can be used here/ They will be added after extensive testing of this program)
        opcWt = getOpctFn(opcFn,opcPts,scalOpc) 
        ### Arrange the Opacity Transfer Function Weights into a Paraview readable 'Points' array
        OPT = opcTFWt(opcTFpts,opcWt,opcPts)
        ### Apply the Opacity weights to the TF (NOTE: This should be a numpy array)
        TFopc.Points = OPT
    else:
        ### Apply the Opacity weights to the TF (NOTE: This should be a numpy array)
        TFopc.Points = opcTbl
    if modColorTF == 1:
        colTF = ColorTF
    ### Processing time calculator - sample 5
    timer5 = time.perf_counter()
    print(f'Calculated Opacity Transfer Function Weights in {round(timer5-timer4, 2)} s')

    ################# Step (OPTIONAL): Apply Text annotation
    if txtAnntn == 1:
        txt_prnt = txt_str + str(format(tmpdat[lp], '.4f'))
        [txt1,txt1D] = addtext(RV1,txt_prnt, txt_ftsz, txt_loc)

    ################# Step 5: Set the camera view
    savFilNam = savrotnam + str(str('%04d' % lp))
    RV1 = set_cam(RV1,rotCam,plnam,ptup,pos_ax,stCampt,rotax,angdiv,savFilNam,hor_res,vert_res,stTranpBG)
 
    ################# Step 6: Save output image
    if rotCam==0:                   # Check if camera rotations were prescribed
        ### Create the file savename
        save_name = savrotnam + str(str('%04d' % lp)) + ".jpg"
        ### Save Screenshot
        SaveScreenshot(save_name, RV1, ImageResolution=[hor_res, vert_res], TransparentBackground=stTranpBG, Quality=100)
        ### Processing time calculator - sample 6
        timer6 = time.perf_counter()
        print(f'Image was saved in {round(timer6-timer5, 2)} s')
    else:                           # If the rotations were applied, they have already been saved
        ### Processing time calculator - sample 6
        timer6 = time.perf_counter()
        print(f'Image Rotations were saved in {round(timer6-timer5, 2)} s') 

    ################# Clear the pipeline
    paraview.simple.Delete(subSetobj)
    del subSetobj
    paraview.simple.Delete(subSetobjDisp)
    del subSetobjDisp
    paraview.simple.Delete(vol)
    del vol
    paraview.simple.Delete(volDisp)
    del volDisp
    paraview.simple.Delete(TF)
    del TF
    paraview.simple.Delete(TFbar)
    del TFbar
    if txtAnntn == 1:
        paraview.simple.Delete(txt1)
        del txt1
        paraview.simple.Delete(txt1D)
        del txt1D
    print('Deleted all objects from Pipeline Browser')
    ResetSession()
    print('Reset the Session')

################# Quit out of the program
exit()
