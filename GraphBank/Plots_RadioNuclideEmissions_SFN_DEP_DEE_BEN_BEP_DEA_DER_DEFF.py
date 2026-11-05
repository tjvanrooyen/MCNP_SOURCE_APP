#####
x_scale = 'linear'
y_scale = 'log'
#####


import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.interpolate import interp1d, PchipInterpolator
from scipy import integrate


# Set default parameters:
plt.rc('lines',  linewidth=7)
plt.rc('figure', figsize=(12,9))
plt.rc('axes',   titlesize=30)
plt.rc('axes',   labelsize=28)
plt.rc('legend', fontsize=24)
plt.rc('xtick',  labelsize=24)
plt.rc('ytick',  labelsize=24)
plt.rc('lines',  markersize=14)

plt.rcParams['figure.constrained_layout.use'] = True
plt.rcParams['figure.constrained_layout.w_pad'] = 0.15   # inches, left+right outer pad (default ~0.042)
plt.rcParams['figure.constrained_layout.h_pad'] = 0.15   # inches, top+bottom outer pad (default ~0.042)
mpl.rcParams['font.family'] = ['Aptos', 'sans-serif']
mpl.rcParams['mathtext.fontset'] = 'dejavuserif'
mpl.rcParams['mathtext.rm'] = 'Aptos'
mpl.rcParams['mathtext.it'] = 'Aptos:italic'
mpl.rcParams['mathtext.bf'] = 'Aptos:bold'


# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "Helvetica"
# })

path_cwd_abs  =  os.getcwd()
path_inp_rel  =  r'Source_Emissions/'
path_out_rel  =  r'GraphBank/'
path_inp_abs  =  os.path.join(path_cwd_abs, path_inp_rel)
path_out_abs  =  os.path.join(path_cwd_abs, path_out_rel)


#####
# Specify consistent line colours:
color_SFN             = 'xkcd:emerald green'
#
color_stemline_DEP    = 'xkcd:true blue'
color_markerface_DEP  = 'xkcd:aqua blue'
color_markeredge_DEP  = 'xkcd:light royal blue'
#
color_stemline_DEE    = 'xkcd:violet'
color_markerface_DEE  = 'xkcd:neon red'
color_markeredge_DEE  = 'xkcd:purply blue'
#
color_BEN             = 'xkcd:blue/purple'
#
color_BEP             = 'xkcd:terracotta'
#
color_stemline_DEA    = 'xkcd:rust brown'
color_markerface_DEA  = 'xkcd:rusty red'
color_markeredge_DEA  = 'xkcd:claret'
#
color_markerface_DER  =  'xkcd:grey'
color_markeredge_DER  =  'xkcd:black'
color_stemline_DER    =  'xkcd:dark grey'

color_grid_major      = 'xkcd:charcoal'
color_grid_minor      = 'xkcd:grey'
color_baseline        = 'xkcd:silver'

color_vlines          = 'xkcd:cherry'
color_vspan           = 'xkcd:pale green'

labelpad_x        = 12
labelpad_y        = 8
position_title_y  = 1.025

#####
lw_grid_major = 1.10
lw_grid_minor = 0.90
lw_baseline = 4
#####

# (1) Identity of RadioNuclide and (2) its Spreadsheet that must be read:


# nuclide = r'Ac223'
# LaTeX_RadioNuclide = r'$\mathrm{^{223}Ac}$'


# nuclide = r'Ac224'
# LaTeX_RadioNuclide = r'$\mathrm{^{224}Ac}$'


# nuclide = r'Ac225'
# LaTeX_RadioNuclide = r'$\mathrm{^{225}Ac}$'


# nuclide = r'Ac225_FE'
# LaTeX_RadioNuclide = r'$\mathrm{^{225}Ac[FE]}$'


# nuclide = r'Ac226'
# LaTeX_RadioNuclide = r'$\mathrm{^{226}Ac}$'


# nuclide = r'Ac227'
# LaTeX_RadioNuclide = r'$\mathrm{^{227}Ac}$'


# nuclide = r'Ac228'
# LaTeX_RadioNuclide = r'$\mathrm{^{228}Ac}$'


# nuclide = r'Ac230'
# LaTeX_RadioNuclide = r'$\mathrm{^{230}Ac}$'


# nuclide = r'Ac231'
# LaTeX_RadioNuclide = r'$\mathrm{^{231}Ac}$'


# nuclide = r'Ac232'
# LaTeX_RadioNuclide = r'$\mathrm{^{232}Ac}$'


# nuclide = r'Ac233'
# LaTeX_RadioNuclide = r'$\mathrm{^{233}Ac}$'


# nuclide = r'Ag99'
# LaTeX_RadioNuclide = r'$\mathrm{^{99}Ag}$'


# nuclide = r'Ag100m'
# LaTeX_RadioNuclide = r'$\mathrm{^{100m}Ag}$'


# nuclide = r'Ag101'
# LaTeX_RadioNuclide = r'$\mathrm{^{101}Ag}$'


# nuclide = r'Ag102'
# LaTeX_RadioNuclide = r'$\mathrm{^{102}Ag}$'


# nuclide = r'Ag102m'
# LaTeX_RadioNuclide = r'$\mathrm{^{102m}Ag}$'


# nuclide = r'Ag103'
# LaTeX_RadioNuclide = r'$\mathrm{^{103}Ag}$'


# nuclide = r'Ag104'
# LaTeX_RadioNuclide = r'$\mathrm{^{104}Ag}$'


# nuclide = r'Ag104m'
# LaTeX_RadioNuclide = r'$\mathrm{^{104m}Ag}$'


# nuclide = r'Ag105'
# LaTeX_RadioNuclide = r'$\mathrm{^{105}Ag}$'


# nuclide = r'Ag105m'
# LaTeX_RadioNuclide = r'$\mathrm{^{105m}Ag}$'


# nuclide = r'Ag106'
# LaTeX_RadioNuclide = r'$\mathrm{^{106}Ag}$'


# nuclide = r'Ag106m'
# LaTeX_RadioNuclide = r'$\mathrm{^{106m}Ag}$'


# nuclide = r'Ag108'
# LaTeX_RadioNuclide = r'$\mathrm{^{108}Ag}$'


# nuclide = r'Ag108m'
# LaTeX_RadioNuclide = r'$\mathrm{^{108m}Ag}$'


# nuclide = r'Ag109m'
# LaTeX_RadioNuclide = r'$\mathrm{^{109m}Ag}$'


# nuclide = r'Ag110'
# LaTeX_RadioNuclide = r'$\mathrm{^{110}Ag}$'


# nuclide = r'Ag110m2'
# LaTeX_RadioNuclide = r'$\mathrm{^{110m2}Ag}$'


# nuclide = r'Ag111'
# LaTeX_RadioNuclide = r'$\mathrm{^{111}Ag}$'


# nuclide = r'Ag111m'
# LaTeX_RadioNuclide = r'$\mathrm{^{111m}Ag}$'


# nuclide = r'Ag112'
# LaTeX_RadioNuclide = r'$\mathrm{^{112}Ag}$'


# nuclide = r'Ag113'
# LaTeX_RadioNuclide = r'$\mathrm{^{113}Ag}$'


# nuclide = r'Ag113m'
# LaTeX_RadioNuclide = r'$\mathrm{^{113m}Ag}$'


# nuclide = r'Ag114'
# LaTeX_RadioNuclide = r'$\mathrm{^{114}Ag}$'


# nuclide = r'Ag115'
# LaTeX_RadioNuclide = r'$\mathrm{^{115}Ag}$'


# nuclide = r'Ag116'
# LaTeX_RadioNuclide = r'$\mathrm{^{116}Ag}$'


# nuclide = r'Ag117'
# LaTeX_RadioNuclide = r'$\mathrm{^{117}Ag}$'


# nuclide = r'Al26'
# LaTeX_RadioNuclide = r'$\mathrm{^{26}Al}$'


# nuclide = r'Al28'
# LaTeX_RadioNuclide = r'$\mathrm{^{28}Al}$'


# nuclide = r'Al29'
# LaTeX_RadioNuclide = r'$\mathrm{^{29}Al}$'


# nuclide = r'Am237'
# LaTeX_RadioNuclide = r'$\mathrm{^{237}Am}$'


# nuclide = r'Am238'
# LaTeX_RadioNuclide = r'$\mathrm{^{238}Am}$'


# nuclide = r'Am239'
# LaTeX_RadioNuclide = r'$\mathrm{^{239}Am}$'


# nuclide = r'Am240'
# LaTeX_RadioNuclide = r'$\mathrm{^{240}Am}$'


# nuclide = r'Am241'
# LaTeX_RadioNuclide = r'$\mathrm{^{241}Am}$'


# nuclide = r'Am241_FE'
# LaTeX_RadioNuclide = r'$\mathrm{^{241}Am[FE]}$'


# nuclide = r'Am242'
# LaTeX_RadioNuclide = r'$\mathrm{^{242}Am}$'


# nuclide = r'Am242m'
# LaTeX_RadioNuclide = r'$\mathrm{^{242m}Am}$'


# nuclide = r'Am243'
# LaTeX_RadioNuclide = r'$\mathrm{^{243}Am}$'


# nuclide = r'Am244'
# LaTeX_RadioNuclide = r'$\mathrm{^{244}Am}$'


# nuclide = r'Am244m'
# LaTeX_RadioNuclide = r'$\mathrm{^{244m}Am}$'


# nuclide = r'Am245'
# LaTeX_RadioNuclide = r'$\mathrm{^{245}Am}$'


# nuclide = r'Am246'
# LaTeX_RadioNuclide = r'$\mathrm{^{246}Am}$'


# nuclide = r'Am246m'
# LaTeX_RadioNuclide = r'$\mathrm{^{246m}Am}$'


# nuclide = r'Am247'
# LaTeX_RadioNuclide = r'$\mathrm{^{247}Am}$'


# nuclide = r'Ar37'
# LaTeX_RadioNuclide = r'$\mathrm{^{37}Ar}$'


# nuclide = r'Ar39'
# LaTeX_RadioNuclide = r'$\mathrm{^{39}Ar}$'


# nuclide = r'Ar41'
# LaTeX_RadioNuclide = r'$\mathrm{^{41}Ar}$'


# nuclide = r'Ar42'
# LaTeX_RadioNuclide = r'$\mathrm{^{42}Ar}$'


# nuclide = r'Ar43'
# LaTeX_RadioNuclide = r'$\mathrm{^{43}Ar}$'


# nuclide = r'Ar44'
# LaTeX_RadioNuclide = r'$\mathrm{^{44}Ar}$'


# nuclide = r'As68'
# LaTeX_RadioNuclide = r'$\mathrm{^{68}As}$'


# nuclide = r'As69'
# LaTeX_RadioNuclide = r'$\mathrm{^{69}As}$'


# nuclide = r'As70'
# LaTeX_RadioNuclide = r'$\mathrm{^{70}As}$'


# nuclide = r'As71'
# LaTeX_RadioNuclide = r'$\mathrm{^{71}As}$'


# nuclide = r'As72'
# LaTeX_RadioNuclide = r'$\mathrm{^{72}As}$'


# nuclide = r'As73'
# LaTeX_RadioNuclide = r'$\mathrm{^{73}As}$'


# nuclide = r'As74'
# LaTeX_RadioNuclide = r'$\mathrm{^{74}As}$'


# nuclide = r'As76'
# LaTeX_RadioNuclide = r'$\mathrm{^{76}As}$'


# nuclide = r'As77'
# LaTeX_RadioNuclide = r'$\mathrm{^{77}As}$'


# nuclide = r'As78'
# LaTeX_RadioNuclide = r'$\mathrm{^{78}As}$'


# nuclide = r'As79'
# LaTeX_RadioNuclide = r'$\mathrm{^{79}As}$'


# nuclide = r'At204'
# LaTeX_RadioNuclide = r'$\mathrm{^{204}At}$'


# nuclide = r'At205'
# LaTeX_RadioNuclide = r'$\mathrm{^{205}At}$'


# nuclide = r'At206'
# LaTeX_RadioNuclide = r'$\mathrm{^{206}At}$'


# nuclide = r'At207'
# LaTeX_RadioNuclide = r'$\mathrm{^{207}At}$'


# nuclide = r'At208'
# LaTeX_RadioNuclide = r'$\mathrm{^{208}At}$'


# nuclide = r'At209'
# LaTeX_RadioNuclide = r'$\mathrm{^{209}At}$'


# nuclide = r'At210'
# LaTeX_RadioNuclide = r'$\mathrm{^{210}At}$'


# nuclide = r'At211'
# LaTeX_RadioNuclide = r'$\mathrm{^{211}At}$'


# nuclide = r'At215'
# LaTeX_RadioNuclide = r'$\mathrm{^{215}At}$'


# nuclide = r'At216'
# LaTeX_RadioNuclide = r'$\mathrm{^{216}At}$'


# nuclide = r'At217'
# LaTeX_RadioNuclide = r'$\mathrm{^{217}At}$'


# nuclide = r'At218'
# LaTeX_RadioNuclide = r'$\mathrm{^{218}At}$'


# nuclide = r'At219'
# LaTeX_RadioNuclide = r'$\mathrm{^{219}At}$'


# nuclide = r'At220'
# LaTeX_RadioNuclide = r'$\mathrm{^{220}At}$'


# nuclide = r'Au186'
# LaTeX_RadioNuclide = r'$\mathrm{^{186}Au}$'


# nuclide = r'Au187'
# LaTeX_RadioNuclide = r'$\mathrm{^{187}Au}$'


# nuclide = r'Au190'
# LaTeX_RadioNuclide = r'$\mathrm{^{190}Au}$'


# nuclide = r'Au191'
# LaTeX_RadioNuclide = r'$\mathrm{^{191}Au}$'


# nuclide = r'Au192'
# LaTeX_RadioNuclide = r'$\mathrm{^{192}Au}$'


# nuclide = r'Au193'
# LaTeX_RadioNuclide = r'$\mathrm{^{193}Au}$'


# nuclide = r'Au193m'
# LaTeX_RadioNuclide = r'$\mathrm{^{193m}Au}$'


# nuclide = r'Au194'
# LaTeX_RadioNuclide = r'$\mathrm{^{194}Au}$'


# nuclide = r'Au195'
# LaTeX_RadioNuclide = r'$\mathrm{^{195}Au}$'


# nuclide = r'Au195m'
# LaTeX_RadioNuclide = r'$\mathrm{^{195m}Au}$'


# nuclide = r'Au196'
# LaTeX_RadioNuclide = r'$\mathrm{^{196}Au}$'


# nuclide = r'Au196m'
# LaTeX_RadioNuclide = r'$\mathrm{^{196m}Au}$'


# nuclide = r'Au198'
# LaTeX_RadioNuclide = r'$\mathrm{^{198}Au}$'


# nuclide = r'Au198m'
# LaTeX_RadioNuclide = r'$\mathrm{^{198m}Au}$'


# nuclide = r'Au199'
# LaTeX_RadioNuclide = r'$\mathrm{^{199}Au}$'


# nuclide = r'Au200'
# LaTeX_RadioNuclide = r'$\mathrm{^{200}Au}$'


# nuclide = r'Au200m'
# LaTeX_RadioNuclide = r'$\mathrm{^{200m}Au}$'


# nuclide = r'Au201'
# LaTeX_RadioNuclide = r'$\mathrm{^{201}Au}$'


# nuclide = r'Au202'
# LaTeX_RadioNuclide = r'$\mathrm{^{202}Au}$'


# nuclide = r'Ba124'
# LaTeX_RadioNuclide = r'$\mathrm{^{124}Ba}$'


# nuclide = r'Ba126'
# LaTeX_RadioNuclide = r'$\mathrm{^{126}Ba}$'


# nuclide = r'Ba127'
# LaTeX_RadioNuclide = r'$\mathrm{^{127}Ba}$'


# nuclide = r'Ba128'
# LaTeX_RadioNuclide = r'$\mathrm{^{128}Ba}$'


# nuclide = r'Ba129'
# LaTeX_RadioNuclide = r'$\mathrm{^{129}Ba}$'


# nuclide = r'Ba129m'
# LaTeX_RadioNuclide = r'$\mathrm{^{129m}Ba}$'


# nuclide = r'Ba131'
# LaTeX_RadioNuclide = r'$\mathrm{^{131}Ba}$'


# nuclide = r'Ba131m'
# LaTeX_RadioNuclide = r'$\mathrm{^{131m}Ba}$'


# nuclide = r'Ba133'
# LaTeX_RadioNuclide = r'$\mathrm{^{133}Ba}$'


# nuclide = r'Ba133m'
# LaTeX_RadioNuclide = r'$\mathrm{^{133m}Ba}$'


# nuclide = r'Ba135m'
# LaTeX_RadioNuclide = r'$\mathrm{^{135m}Ba}$'


# nuclide = r'Ba137m'
# LaTeX_RadioNuclide = r'$\mathrm{^{137m}Ba}$'


# nuclide = r'Ba139'
# LaTeX_RadioNuclide = r'$\mathrm{^{139}Ba}$'


# nuclide = r'Ba140'
# LaTeX_RadioNuclide = r'$\mathrm{^{140}Ba}$'


# nuclide = r'Ba141'
# LaTeX_RadioNuclide = r'$\mathrm{^{141}Ba}$'


# nuclide = r'Ba142'
# LaTeX_RadioNuclide = r'$\mathrm{^{142}Ba}$'


# nuclide = r'Be7'
# LaTeX_RadioNuclide = r'$\mathrm{^{7}Be}$'


# nuclide = r'Be10'
# LaTeX_RadioNuclide = r'$\mathrm{^{10}Be}$'


# nuclide = r'Bi197'
# LaTeX_RadioNuclide = r'$\mathrm{^{197}Bi}$'


# nuclide = r'Bi200'
# LaTeX_RadioNuclide = r'$\mathrm{^{200}Bi}$'


# nuclide = r'Bi201'
# LaTeX_RadioNuclide = r'$\mathrm{^{201}Bi}$'


# nuclide = r'Bi202'
# LaTeX_RadioNuclide = r'$\mathrm{^{202}Bi}$'


# nuclide = r'Bi203'
# LaTeX_RadioNuclide = r'$\mathrm{^{203}Bi}$'


# nuclide = r'Bi204'
# LaTeX_RadioNuclide = r'$\mathrm{^{204}Bi}$'


# nuclide = r'Bi205'
# LaTeX_RadioNuclide = r'$\mathrm{^{205}Bi}$'


# nuclide = r'Bi206'
# LaTeX_RadioNuclide = r'$\mathrm{^{206}Bi}$'


# nuclide = r'Bi207'
# LaTeX_RadioNuclide = r'$\mathrm{^{207}Bi}$'


# nuclide = r'Bi208'
# LaTeX_RadioNuclide = r'$\mathrm{^{208}Bi}$'


# nuclide = r'Bi210'
# LaTeX_RadioNuclide = r'$\mathrm{^{210}Bi}$'


# nuclide = r'Bi210m'
# LaTeX_RadioNuclide = r'$\mathrm{^{210m}Bi}$'


# nuclide = r'Bi211'
# LaTeX_RadioNuclide = r'$\mathrm{^{211}Bi}$'


# nuclide = r'Bi212'
# LaTeX_RadioNuclide = r'$\mathrm{^{212}Bi}$'


# nuclide = r'Bi212m2'
# LaTeX_RadioNuclide = r'$\mathrm{^{212m2}Bi}$'


# nuclide = r'Bi213'
# LaTeX_RadioNuclide = r'$\mathrm{^{213}Bi}$'


# nuclide = r'Bi214'
# LaTeX_RadioNuclide = r'$\mathrm{^{214}Bi}$'


# nuclide = r'Bi215'
# LaTeX_RadioNuclide = r'$\mathrm{^{215}Bi}$'


# nuclide = r'Bi216'
# LaTeX_RadioNuclide = r'$\mathrm{^{216}Bi}$'


# nuclide = r'Bk245'
# LaTeX_RadioNuclide = r'$\mathrm{^{245}Bk}$'


# nuclide = r'Bk246'
# LaTeX_RadioNuclide = r'$\mathrm{^{246}Bk}$'


# nuclide = r'Bk247'
# LaTeX_RadioNuclide = r'$\mathrm{^{247}Bk}$'


# nuclide = r'Bk248m'
# LaTeX_RadioNuclide = r'$\mathrm{^{248m}Bk}$'


# nuclide = r'Bk249'
# LaTeX_RadioNuclide = r'$\mathrm{^{249}Bk}$'


# nuclide = r'Bk250'
# LaTeX_RadioNuclide = r'$\mathrm{^{250}Bk}$'


# nuclide = r'Bk251'
# LaTeX_RadioNuclide = r'$\mathrm{^{251}Bk}$'


# nuclide = r'Br72'
# LaTeX_RadioNuclide = r'$\mathrm{^{72}Br}$'


# nuclide = r'Br73'
# LaTeX_RadioNuclide = r'$\mathrm{^{73}Br}$'


# nuclide = r'Br74'
# LaTeX_RadioNuclide = r'$\mathrm{^{74}Br}$'


# nuclide = r'Br74m'
# LaTeX_RadioNuclide = r'$\mathrm{^{74m}Br}$'


# nuclide = r'Br75'
# LaTeX_RadioNuclide = r'$\mathrm{^{75}Br}$'


# nuclide = r'Br76'
# LaTeX_RadioNuclide = r'$\mathrm{^{76}Br}$'


# nuclide = r'Br76m'
# LaTeX_RadioNuclide = r'$\mathrm{^{76m}Br}$'


# nuclide = r'Br77'
# LaTeX_RadioNuclide = r'$\mathrm{^{77}Br}$'


# nuclide = r'Br77m'
# LaTeX_RadioNuclide = r'$\mathrm{^{77m}Br}$'


# nuclide = r'Br78'
# LaTeX_RadioNuclide = r'$\mathrm{^{78}Br}$'


# nuclide = r'Br80'
# LaTeX_RadioNuclide = r'$\mathrm{^{80}Br}$'


# nuclide = r'Br80m'
# LaTeX_RadioNuclide = r'$\mathrm{^{80m}Br}$'


# nuclide = r'Br82'
# LaTeX_RadioNuclide = r'$\mathrm{^{82}Br}$'


# nuclide = r'Br82m'
# LaTeX_RadioNuclide = r'$\mathrm{^{82m}Br}$'


# nuclide = r'Br83'
# LaTeX_RadioNuclide = r'$\mathrm{^{83}Br}$'


# nuclide = r'Br84'
# LaTeX_RadioNuclide = r'$\mathrm{^{84}Br}$'


# nuclide = r'Br84m'
# LaTeX_RadioNuclide = r'$\mathrm{^{84m}Br}$'


# nuclide = r'Br85'
# LaTeX_RadioNuclide = r'$\mathrm{^{85}Br}$'


# nuclide = r'C10'
# LaTeX_RadioNuclide = r'$\mathrm{^{10}C}$'


# nuclide = r'C11'
# LaTeX_RadioNuclide = r'$\mathrm{^{11}C}$'


# nuclide = r'C14'
# LaTeX_RadioNuclide = r'$\mathrm{^{14}C}$'


# nuclide = r'Ca41'
# LaTeX_RadioNuclide = r'$\mathrm{^{41}Ca}$'


# nuclide = r'Ca45'
# LaTeX_RadioNuclide = r'$\mathrm{^{45}Ca}$'


# nuclide = r'Ca47'
# LaTeX_RadioNuclide = r'$\mathrm{^{47}Ca}$'


# nuclide = r'Ca49'
# LaTeX_RadioNuclide = r'$\mathrm{^{49}Ca}$'


# nuclide = r'Cd101'
# LaTeX_RadioNuclide = r'$\mathrm{^{101}Cd}$'


# nuclide = r'Cd102'
# LaTeX_RadioNuclide = r'$\mathrm{^{102}Cd}$'


# nuclide = r'Cd103'
# LaTeX_RadioNuclide = r'$\mathrm{^{103}Cd}$'


# nuclide = r'Cd104'
# LaTeX_RadioNuclide = r'$\mathrm{^{104}Cd}$'


# nuclide = r'Cd105'
# LaTeX_RadioNuclide = r'$\mathrm{^{105}Cd}$'


# nuclide = r'Cd107'
# LaTeX_RadioNuclide = r'$\mathrm{^{107}Cd}$'


# nuclide = r'Cd109'
# LaTeX_RadioNuclide = r'$\mathrm{^{109}Cd}$'


# nuclide = r'Cd111m'
# LaTeX_RadioNuclide = r'$\mathrm{^{111m}Cd}$'


# nuclide = r'Cd113'
# LaTeX_RadioNuclide = r'$\mathrm{^{113}Cd}$'


# nuclide = r'Cd113m'
# LaTeX_RadioNuclide = r'$\mathrm{^{113m}Cd}$'


# nuclide = r'Cd115'
# LaTeX_RadioNuclide = r'$\mathrm{^{115}Cd}$'


# nuclide = r'Cd115m'
# LaTeX_RadioNuclide = r'$\mathrm{^{115m}Cd}$'


# nuclide = r'Cd117'
# LaTeX_RadioNuclide = r'$\mathrm{^{117}Cd}$'


# nuclide = r'Cd117m'
# LaTeX_RadioNuclide = r'$\mathrm{^{117m}Cd}$'


# nuclide = r'Cd118'
# LaTeX_RadioNuclide = r'$\mathrm{^{118}Cd}$'


# nuclide = r'Cd119'
# LaTeX_RadioNuclide = r'$\mathrm{^{119}Cd}$'


# nuclide = r'Cd119m'
# LaTeX_RadioNuclide = r'$\mathrm{^{119m}Cd}$'


# nuclide = r'Ce130'
# LaTeX_RadioNuclide = r'$\mathrm{^{130}Ce}$'


# nuclide = r'Ce131'
# LaTeX_RadioNuclide = r'$\mathrm{^{131}Ce}$'


# nuclide = r'Ce132'
# LaTeX_RadioNuclide = r'$\mathrm{^{132}Ce}$'


# nuclide = r'Ce133'
# LaTeX_RadioNuclide = r'$\mathrm{^{133}Ce}$'


# nuclide = r'Ce133m'
# LaTeX_RadioNuclide = r'$\mathrm{^{133m}Ce}$'


# nuclide = r'Ce134'
# LaTeX_RadioNuclide = r'$\mathrm{^{134}Ce}$'


# nuclide = r'Ce135'
# LaTeX_RadioNuclide = r'$\mathrm{^{135}Ce}$'


# nuclide = r'Ce137'
# LaTeX_RadioNuclide = r'$\mathrm{^{137}Ce}$'


# nuclide = r'Ce137m'
# LaTeX_RadioNuclide = r'$\mathrm{^{137m}Ce}$'


# nuclide = r'Ce139'
# LaTeX_RadioNuclide = r'$\mathrm{^{139}Ce}$'


# nuclide = r'Ce141'
# LaTeX_RadioNuclide = r'$\mathrm{^{141}Ce}$'


# nuclide = r'Ce143'
# LaTeX_RadioNuclide = r'$\mathrm{^{143}Ce}$'


# nuclide = r'Ce144'
# LaTeX_RadioNuclide = r'$\mathrm{^{144}Ce}$'


# nuclide = r'Ce144_Pr144'
# LaTeX_RadioNuclide = r'$\mathrm{^{144m}Ce + ^{144}Pr}$'


# nuclide = r'Ce145'
# LaTeX_RadioNuclide = r'$\mathrm{^{145}Ce}$'


# nuclide = r'Cf244'
# LaTeX_RadioNuclide = r'$\mathrm{^{244}Cf}$'


# nuclide = r'Cf246'
# LaTeX_RadioNuclide = r'$\mathrm{^{246}Cf}$'


# nuclide = r'Cf247'
# LaTeX_RadioNuclide = r'$\mathrm{^{247}Cf}$'


# nuclide = r'Cf248'
# LaTeX_RadioNuclide = r'$\mathrm{^{248}Cf}$'


# nuclide = r'Cf249'
# LaTeX_RadioNuclide = r'$\mathrm{^{249}Cf}$'


# nuclide = r'Cf250'
# LaTeX_RadioNuclide = r'$\mathrm{^{250}Cf}$'


# nuclide = r'Cf251'
# LaTeX_RadioNuclide = r'$\mathrm{^{251}Cf}$'


# nuclide = r'Cf252'
# LaTeX_RadioNuclide = r'$\mathrm{^{252}Cf}$'


# nuclide = r'Cf252_FE'
# LaTeX_RadioNuclide = r'$\mathrm{^{252}Cf[FE]}$'


# nuclide = r'Cf253'
# LaTeX_RadioNuclide = r'$\mathrm{^{253}Cf}$'


# nuclide = r'Cf254'
# LaTeX_RadioNuclide = r'$\mathrm{^{254}Cf}$'


# nuclide = r'Cf255'
# LaTeX_RadioNuclide = r'$\mathrm{^{255}Cf}$'


# nuclide = r'Cl34'
# LaTeX_RadioNuclide = r'$\mathrm{^{34}Cl}$'


# nuclide = r'Cl34m'
# LaTeX_RadioNuclide = r'$\mathrm{^{34m}Cl}$'


# nuclide = r'Cl36'
# LaTeX_RadioNuclide = r'$\mathrm{_{17}^{36}Cl}$'


# nuclide = r'Cl38'
# LaTeX_RadioNuclide = r'$\mathrm{^{38}Cl}$'


# nuclide = r'Cl39'
# LaTeX_RadioNuclide = r'$\mathrm{^{39}Cl}$'


# nuclide = r'Cl40'
# LaTeX_RadioNuclide = r'$\mathrm{^{40}Cl}$'


# nuclide = r'Cm238'
# LaTeX_RadioNuclide = r'$\mathrm{^{238}Cm}$'


# nuclide = r'Cm239'
# LaTeX_RadioNuclide = r'$\mathrm{^{239}Cm}$'


# nuclide = r'Cm240'
# LaTeX_RadioNuclide = r'$\mathrm{^{240}Cm}$'


# nuclide = r'Cm241'
# LaTeX_RadioNuclide = r'$\mathrm{^{241}Cm}$'


# nuclide = r'Cm242'
# LaTeX_RadioNuclide = r'$\mathrm{^{242}Cm}$'


# nuclide = r'Cm243'
# LaTeX_RadioNuclide = r'$\mathrm{^{243}Cm}$'


# nuclide = r'Cm244'
# LaTeX_RadioNuclide = r'$\mathrm{^{244}Cm}$'


# nuclide = r'Cm245'
# LaTeX_RadioNuclide = r'$\mathrm{^{245}Cm}$'


# nuclide = r'Cm246'
# LaTeX_RadioNuclide = r'$\mathrm{^{246}Cm}$'


# nuclide = r'Cm247'
# LaTeX_RadioNuclide = r'$\mathrm{^{247}Cm}$'


# nuclide = r'Cm248'
# LaTeX_RadioNuclide = r'$\mathrm{^{248}Cm}$'


# nuclide = r'Cm249'
# LaTeX_RadioNuclide = r'$\mathrm{^{249}Cm}$'


# nuclide = r'Cm250'
# LaTeX_RadioNuclide = r'$\mathrm{^{250}Cm}$'


# nuclide = r'Cm251'
# LaTeX_RadioNuclide = r'$\mathrm{^{251}Cm}$'


# nuclide = r'Co54m'
# LaTeX_RadioNuclide = r'$\mathrm{^{54m}Co}$'


# nuclide = r'Co55'
# LaTeX_RadioNuclide = r'$\mathrm{^{55}Co}$'


# nuclide = r'Co56'
# LaTeX_RadioNuclide = r'$\mathrm{_{27}^{56}Co}$'


# nuclide = r'Co57'
# LaTeX_RadioNuclide = r'$\mathrm{_{27}^{57}Co}$'


# nuclide = r'Co58'
# LaTeX_RadioNuclide = r'$\mathrm{_{27}^{58}Co}$'


# nuclide = r'Co58m'
# LaTeX_RadioNuclide = r'$\mathrm{^{58m}Co}$'


# nuclide = r'Co60'
# LaTeX_RadioNuclide = r'$\mathrm{_{27}^{60}Co}$'


# nuclide = r'Co60m'
# LaTeX_RadioNuclide = r'$\mathrm{^{60m}Co}$'


# nuclide = r'Co61'
# LaTeX_RadioNuclide = r'$\mathrm{^{61}Co}$'


# nuclide = r'Co62'
# LaTeX_RadioNuclide = r'$\mathrm{^{62}Co}$'


# nuclide = r'Co62m'
# LaTeX_RadioNuclide = r'$\mathrm{^{62m}Co}$'


# nuclide = r'Cr48'
# LaTeX_RadioNuclide = r'$\mathrm{^{48}Cr}$'


# nuclide = r'Cr49'
# LaTeX_RadioNuclide = r'$\mathrm{^{49}Cr}$'


# nuclide = r'Cr51'
# LaTeX_RadioNuclide = r'$\mathrm{_{24}^{51}Cr}$'


# nuclide = r'Cr55'
# LaTeX_RadioNuclide = r'$\mathrm{^{55}Cr}$'


# nuclide = r'Cr56'
# LaTeX_RadioNuclide = r'$\mathrm{^{56}Cr}$'


# nuclide = r'Cs121'
# LaTeX_RadioNuclide = r'$\mathrm{^{121}Cs}$'


# nuclide = r'Cs121m'
# LaTeX_RadioNuclide = r'$\mathrm{^{121m}Cs}$'


# nuclide = r'Cs123'
# LaTeX_RadioNuclide = r'$\mathrm{^{123}Cs}$'


# nuclide = r'Cs124'
# LaTeX_RadioNuclide = r'$\mathrm{^{124}Cs}$'


# nuclide = r'Cs125'
# LaTeX_RadioNuclide = r'$\mathrm{^{125}Cs}$'


# nuclide = r'Cs126'
# LaTeX_RadioNuclide = r'$\mathrm{^{126}Cs}$'


# nuclide = r'Cs127'
# LaTeX_RadioNuclide = r'$\mathrm{^{127}Cs}$'


# nuclide = r'Cs128'
# LaTeX_RadioNuclide = r'$\mathrm{^{128}Cs}$'


# nuclide = r'Cs129'
# LaTeX_RadioNuclide = r'$\mathrm{^{129}Cs}$'


# nuclide = r'Cs130'
# LaTeX_RadioNuclide = r'$\mathrm{^{130}Cs}$'


# nuclide = r'Cs130m'
# LaTeX_RadioNuclide = r'$\mathrm{^{130m}Cs}$'


# nuclide = r'Cs131'
# LaTeX_RadioNuclide = r'$\mathrm{^{131}Cs}$'


# nuclide = r'Cs132'
# LaTeX_RadioNuclide = r'$\mathrm{^{132}Cs}$'


# nuclide = r'Cs134'
# LaTeX_RadioNuclide = r'$\mathrm{^{134}Cs}$'


# nuclide = r'Cs134m'
# LaTeX_RadioNuclide = r'$\mathrm{^{134m}Cs}$'


# nuclide = r'Cs135'
# LaTeX_RadioNuclide = r'$\mathrm{^{135}Cs}$'


# nuclide = r'Cs135m'
# LaTeX_RadioNuclide = r'$\mathrm{^{135m}Cs}$'


# nuclide = r'Cs136'
# LaTeX_RadioNuclide = r'$\mathrm{^{136}Cs}$'


# nuclide = r'Cs137'
# LaTeX_RadioNuclide = r'$\mathrm{^{137}Cs}$'


# nuclide = r'Cs137_Ba137m'
# LaTeX_RadioNuclide = r'$\mathrm{^{137}Cs} + \mathrm{^{137m}Ba}$'


# nuclide = r'Cs138'
# LaTeX_RadioNuclide = r'$\mathrm{^{138}Cs}$'


# nuclide = r'Cs138m'
# LaTeX_RadioNuclide = r'$\mathrm{^{138m}Cs}$'


# nuclide = r'Cs139'
# LaTeX_RadioNuclide = r'$\mathrm{^{139}Cs}$'


# nuclide = r'Cs140'
# LaTeX_RadioNuclide = r'$\mathrm{^{140}Cs}$'


# nuclide = r'Cu57'
# LaTeX_RadioNuclide = r'$\mathrm{^{57}Cu}$'


# nuclide = r'Cu59'
# LaTeX_RadioNuclide = r'$\mathrm{^{59}Cu}$'


# nuclide = r'Cu60'
# LaTeX_RadioNuclide = r'$\mathrm{^{60}Cu}$'


# nuclide = r'Cu61'
# LaTeX_RadioNuclide = r'$\mathrm{^{61}Cu}$'


# nuclide = r'Cu62'
# LaTeX_RadioNuclide = r'$\mathrm{^{62}Cu}$'


# nuclide = r'Cu64'
# LaTeX_RadioNuclide = r'$\mathrm{_{29}^{64}Cu}$'


# nuclide = r'Cu66'
# LaTeX_RadioNuclide = r'$\mathrm{^{66}Cu}$'


# nuclide = r'Cu67'
# LaTeX_RadioNuclide = r'$\mathrm{_{29}^{67}Cu}$'


# nuclide = r'Cu69'
# LaTeX_RadioNuclide = r'$\mathrm{^{69}Cu}$'


# nuclide = r'Dy148'
# LaTeX_RadioNuclide = r'$\mathrm{^{148}Dy}$'


# nuclide = r'Dy149'
# LaTeX_RadioNuclide = r'$\mathrm{^{149}Dy}$'


# nuclide = r'Dy150'
# LaTeX_RadioNuclide = r'$\mathrm{^{150}Dy}$'


# nuclide = r'Dy151'
# LaTeX_RadioNuclide = r'$\mathrm{^{151}Dy}$'


# nuclide = r'Dy152'
# LaTeX_RadioNuclide = r'$\mathrm{^{152}Dy}$'


# nuclide = r'Dy153'
# LaTeX_RadioNuclide = r'$\mathrm{^{153}Dy}$'


# nuclide = r'Dy154'
# LaTeX_RadioNuclide = r'$\mathrm{^{154}Dy}$'


# nuclide = r'Dy155'
# LaTeX_RadioNuclide = r'$\mathrm{^{155}Dy}$'


# nuclide = r'Dy157'
# LaTeX_RadioNuclide = r'$\mathrm{^{157}Dy}$'


# nuclide = r'Dy159'
# LaTeX_RadioNuclide = r'$\mathrm{^{159}Dy}$'


# nuclide = r'Dy165'
# LaTeX_RadioNuclide = r'$\mathrm{^{165}Dy}$'


# nuclide = r'Dy165m'
# LaTeX_RadioNuclide = r'$\mathrm{^{165m}Dy}$'


# nuclide = r'Dy166'
# LaTeX_RadioNuclide = r'$\mathrm{^{166}Dy}$'


# nuclide = r'Dy167'
# LaTeX_RadioNuclide = r'$\mathrm{^{167}Dy}$'


# nuclide = r'Dy168'
# LaTeX_RadioNuclide = r'$\mathrm{^{168}Dy}$'


# nuclide = r'Er154'
# LaTeX_RadioNuclide = r'$\mathrm{^{154}Er}$'


# nuclide = r'Er156'
# LaTeX_RadioNuclide = r'$\mathrm{^{156}Er}$'


# nuclide = r'Er159'
# LaTeX_RadioNuclide = r'$\mathrm{^{159}Er}$'


# nuclide = r'Er161'
# LaTeX_RadioNuclide = r'$\mathrm{^{161}Er}$'


# nuclide = r'Er163'
# LaTeX_RadioNuclide = r'$\mathrm{^{163}Er}$'


# nuclide = r'Er165'
# LaTeX_RadioNuclide = r'$\mathrm{^{165}Er}$'


# nuclide = r'Er167m'
# LaTeX_RadioNuclide = r'$\mathrm{^{167m}Er}$'


# nuclide = r'Er169'
# LaTeX_RadioNuclide = r'$\mathrm{^{169}Er}$'


# nuclide = r'Er171'
# LaTeX_RadioNuclide = r'$\mathrm{^{171}Er}$'


# nuclide = r'Er172'
# LaTeX_RadioNuclide = r'$\mathrm{^{172}Er}$'


# nuclide = r'Er173'
# LaTeX_RadioNuclide = r'$\mathrm{^{173}Er}$'


# nuclide = r'Es249'
# LaTeX_RadioNuclide = r'$\mathrm{^{249}Es}$'


# nuclide = r'Es250'
# LaTeX_RadioNuclide = r'$\mathrm{^{250}Es}$'


# nuclide = r'Es250m'
# LaTeX_RadioNuclide = r'$\mathrm{^{250m}Es}$'


# nuclide = r'Es251'
# LaTeX_RadioNuclide = r'$\mathrm{^{251}Es}$'


# nuclide = r'Es253'
# LaTeX_RadioNuclide = r'$\mathrm{^{253}Es}$'


# nuclide = r'Es254'
# LaTeX_RadioNuclide = r'$\mathrm{^{254}Es}$'


# nuclide = r'Es254m'
# LaTeX_RadioNuclide = r'$\mathrm{^{254m}Es}$'


# nuclide = r'Es255'
# LaTeX_RadioNuclide = r'$\mathrm{^{255}Es}$'


# nuclide = r'Es256'
# LaTeX_RadioNuclide = r'$\mathrm{^{256}Es}$'


# nuclide = r'Eu142'
# LaTeX_RadioNuclide = r'$\mathrm{^{142}Eu}$'


# nuclide = r'Eu142m'
# LaTeX_RadioNuclide = r'$\mathrm{^{142m}Eu}$'


# nuclide = r'Eu143'
# LaTeX_RadioNuclide = r'$\mathrm{^{143}Eu}$'


# nuclide = r'Eu144'
# LaTeX_RadioNuclide = r'$\mathrm{^{144}Eu}$'


# nuclide = r'Eu145'
# LaTeX_RadioNuclide = r'$\mathrm{^{145}Eu}$'


# nuclide = r'Eu146'
# LaTeX_RadioNuclide = r'$\mathrm{^{146}Eu}$'


# nuclide = r'Eu147'
# LaTeX_RadioNuclide = r'$\mathrm{^{147}Eu}$'


# nuclide = r'Eu148'
# LaTeX_RadioNuclide = r'$\mathrm{^{148}Eu}$'


# nuclide = r'Eu149'
# LaTeX_RadioNuclide = r'$\mathrm{^{149}Eu}$'


# nuclide = r'Eu150'
# LaTeX_RadioNuclide = r'$\mathrm{^{150}Eu}$'


# nuclide = r'Eu150m'
# LaTeX_RadioNuclide = r'$\mathrm{^{150m}Eu}$'


# nuclide = r'Eu152'
# LaTeX_RadioNuclide = r'$\mathrm{^{152}Eu}$'


# nuclide = r'Eu152m'
# LaTeX_RadioNuclide = r'$\mathrm{^{152m}Eu}$'


# nuclide = r'Eu152m2'
# LaTeX_RadioNuclide = r'$\mathrm{^{152m2}Eu}$'


# nuclide = r'Eu154'
# LaTeX_RadioNuclide = r'$\mathrm{^{154}Eu}$'


# nuclide = r'Eu154m'
# LaTeX_RadioNuclide = r'$\mathrm{^{154m}Eu}$'


# nuclide = r'Eu155'
# LaTeX_RadioNuclide = r'$\mathrm{^{155}Eu}$'


# nuclide = r'Eu156'
# LaTeX_RadioNuclide = r'$\mathrm{^{156}Eu}$'


# nuclide = r'Eu157'
# LaTeX_RadioNuclide = r'$\mathrm{^{157}Eu}$'


# nuclide = r'Eu158'
# LaTeX_RadioNuclide = r'$\mathrm{^{158}Eu}$'


# nuclide = r'Eu159'
# LaTeX_RadioNuclide = r'$\mathrm{^{159}Eu}$'


# nuclide = r'F17'
# LaTeX_RadioNuclide = r'$\mathrm{^{17}F}$'


# nuclide = r'F18'
# LaTeX_RadioNuclide = r'$\mathrm{^{18}F}$'


# nuclide = r'Fe52'
# LaTeX_RadioNuclide = r'$\mathrm{^{52}Fe}$'


# nuclide = r'Fe53'
# LaTeX_RadioNuclide = r'$\mathrm{^{53}Fe}$'


# nuclide = r'Fe53m'
# LaTeX_RadioNuclide = r'$\mathrm{^{53m}Fe}$'


# nuclide = r'Fe55'
# LaTeX_RadioNuclide = r'$\mathrm{_{26}^{55}Fe}$'


# nuclide = r'Fe59'
# LaTeX_RadioNuclide = r'$\mathrm{_{26}^{59}Fe}$'


# nuclide = r'Fe60'
# LaTeX_RadioNuclide = r'$\mathrm{^{60}Fe}$'


# nuclide = r'Fe61'
# LaTeX_RadioNuclide = r'$\mathrm{^{61}Fe}$'


# nuclide = r'Fe62'
# LaTeX_RadioNuclide = r'$\mathrm{^{62}Fe}$'


# nuclide = r'Fm251'
# LaTeX_RadioNuclide = r'$\mathrm{^{251}Fm}$'


# nuclide = r'Fm252'
# LaTeX_RadioNuclide = r'$\mathrm{^{252}Fm}$'


# nuclide = r'Fm253'
# LaTeX_RadioNuclide = r'$\mathrm{^{253}Fm}$'


# nuclide = r'Fm254'
# LaTeX_RadioNuclide = r'$\mathrm{^{254}Fm}$'


# nuclide = r'Fm255'
# LaTeX_RadioNuclide = r'$\mathrm{^{255}Fm}$'


# nuclide = r'Fm256'
# LaTeX_RadioNuclide = r'$\mathrm{^{256}Fm}$'


# nuclide = r'Fm257'
# LaTeX_RadioNuclide = r'$\mathrm{^{257}Fm}$'


# nuclide = r'Fr212'
# LaTeX_RadioNuclide = r'$\mathrm{^{212}Fr}$'


# nuclide = r'Fr219'
# LaTeX_RadioNuclide = r'$\mathrm{^{219}Fr}$'


# nuclide = r'Fr220'
# LaTeX_RadioNuclide = r'$\mathrm{^{220}Fr}$'


# nuclide = r'Fr221'
# LaTeX_RadioNuclide = r'$\mathrm{^{221}Fr}$'


# nuclide = r'Fr222'
# LaTeX_RadioNuclide = r'$\mathrm{^{222}Fr}$'


# nuclide = r'Fr223'
# LaTeX_RadioNuclide = r'$\mathrm{^{223}Fr}$'


# nuclide = r'Fr224'
# LaTeX_RadioNuclide = r'$\mathrm{^{224}Fr}$'


# nuclide = r'Fr227'
# LaTeX_RadioNuclide = r'$\mathrm{^{227}Fr}$'


# nuclide = r'Ga64'
# LaTeX_RadioNuclide = r'$\mathrm{^{64}Ga}$'


# nuclide = r'Ga65'
# LaTeX_RadioNuclide = r'$\mathrm{^{65}Ga}$'


# nuclide = r'Ga66'
# LaTeX_RadioNuclide = r'$\mathrm{^{66}Ga}$'


# nuclide = r'Ga67'
# LaTeX_RadioNuclide = r'$\mathrm{_{31}^{67}Ga}$'


# nuclide = r'Ga68'
# LaTeX_RadioNuclide = r'$\mathrm{_{31}^{68}Ga}$'


# nuclide = r'Ga70'
# LaTeX_RadioNuclide = r'$\mathrm{^{70}Ga}$'


# nuclide = r'Ga72'
# LaTeX_RadioNuclide = r'$\mathrm{^{72}Ga}$'


# nuclide = r'Ga73'
# LaTeX_RadioNuclide = r'$\mathrm{^{73}Ga}$'


# nuclide = r'Ga74'
# LaTeX_RadioNuclide = r'$\mathrm{^{74}Ga}$'


# nuclide = r'Gd142'
# LaTeX_RadioNuclide = r'$\mathrm{^{142}Gd}$'


# nuclide = r'Gd143m'
# LaTeX_RadioNuclide = r'$\mathrm{^{143m}Gd}$'


# nuclide = r'Gd144'
# LaTeX_RadioNuclide = r'$\mathrm{^{144}Gd}$'


# nuclide = r'Gd145'
# LaTeX_RadioNuclide = r'$\mathrm{^{145}Gd}$'


# nuclide = r'Gd145m'
# LaTeX_RadioNuclide = r'$\mathrm{^{145m}Gd}$'


# nuclide = r'Gd146'
# LaTeX_RadioNuclide = r'$\mathrm{^{146}Gd}$'


# nuclide = r'Gd147'
# LaTeX_RadioNuclide = r'$\mathrm{^{147}Gd}$'


# nuclide = r'Gd148'
# LaTeX_RadioNuclide = r'$\mathrm{^{148}Gd}$'


# nuclide = r'Gd149'
# LaTeX_RadioNuclide = r'$\mathrm{^{149}Gd}$'


# nuclide = r'Gd150'
# LaTeX_RadioNuclide = r'$\mathrm{^{150}Gd}$'


# nuclide = r'Gd151'
# LaTeX_RadioNuclide = r'$\mathrm{^{151}Gd}$'


# nuclide = r'Gd152'
# LaTeX_RadioNuclide = r'$\mathrm{^{152}Gd}$'


# nuclide = r'Gd153'
# LaTeX_RadioNuclide = r'$\mathrm{^{153}Gd}$'


# nuclide = r'Gd159'
# LaTeX_RadioNuclide = r'$\mathrm{^{159}Gd}$'


# nuclide = r'Gd162'
# LaTeX_RadioNuclide = r'$\mathrm{^{162}Gd}$'


# nuclide = r'Ge66'
# LaTeX_RadioNuclide = r'$\mathrm{^{66}Ge}$'


# nuclide = r'Ge67'
# LaTeX_RadioNuclide = r'$\mathrm{^{67}Ge}$'


# nuclide = r'Ge68'
# LaTeX_RadioNuclide = r'$\mathrm{_{32}^{68}Ge}$'


# nuclide = r'Ge69'
# LaTeX_RadioNuclide = r'$\mathrm{^{69}Ge}$'


# nuclide = r'Ge71'
# LaTeX_RadioNuclide = r'$\mathrm{^{71}Ge}$'


# nuclide = r'Ge75'
# LaTeX_RadioNuclide = r'$\mathrm{^{75}Ge}$'


# nuclide = r'Ge77'
# LaTeX_RadioNuclide = r'$\mathrm{^{77}Ge}$'


# nuclide = r'Ge78'
# LaTeX_RadioNuclide = r'$\mathrm{^{78}Ge}$'


# nuclide = r'H3'
# LaTeX_RadioNuclide = r'$\mathrm{_{1}^{3}H}$'


# nuclide = r'Hf167'
# LaTeX_RadioNuclide = r'$\mathrm{^{167}Hf}$'


# nuclide = r'Hf169'
# LaTeX_RadioNuclide = r'$\mathrm{^{169}Hf}$'


# nuclide = r'Hf170'
# LaTeX_RadioNuclide = r'$\mathrm{^{170}Hf}$'


# nuclide = r'Hf172'
# LaTeX_RadioNuclide = r'$\mathrm{^{172}Hf}$'


# nuclide = r'Hf173'
# LaTeX_RadioNuclide = r'$\mathrm{^{173}Hf}$'


# nuclide = r'Hf174'
# LaTeX_RadioNuclide = r'$\mathrm{^{174}Hf}$'


# nuclide = r'Hf175'
# LaTeX_RadioNuclide = r'$\mathrm{^{175}Hf}$'


# nuclide = r'Hf177m'
# LaTeX_RadioNuclide = r'$\mathrm{^{177m}Hf}$'


# nuclide = r'Hf178m'
# LaTeX_RadioNuclide = r'$\mathrm{^{178m}Hf}$'


# nuclide = r'Hf179m'
# LaTeX_RadioNuclide = r'$\mathrm{^{179m}Hf}$'


# nuclide = r'Hf180m'
# LaTeX_RadioNuclide = r'$\mathrm{^{180m}Hf}$'


# nuclide = r'Hf181'
# LaTeX_RadioNuclide = r'$\mathrm{^{181}Hf}$'


# nuclide = r'Hf182'
# LaTeX_RadioNuclide = r'$\mathrm{^{182}Hf}$'


# nuclide = r'Hf182m'
# LaTeX_RadioNuclide = r'$\mathrm{^{182m}Hf}$'


# nuclide = r'Hf183'
# LaTeX_RadioNuclide = r'$\mathrm{^{183}Hf}$'


# nuclide = r'Hf184'
# LaTeX_RadioNuclide = r'$\mathrm{^{184}Hf}$'


# nuclide = r'Hg190'
# LaTeX_RadioNuclide = r'$\mathrm{^{190}Hg}$'


# nuclide = r'Hg191m'
# LaTeX_RadioNuclide = r'$\mathrm{^{191m}Hg}$'


# nuclide = r'Hg192'
# LaTeX_RadioNuclide = r'$\mathrm{^{192}Hg}$'


# nuclide = r'Hg193'
# LaTeX_RadioNuclide = r'$\mathrm{^{193}Hg}$'


# nuclide = r'Hg193m'
# LaTeX_RadioNuclide = r'$\mathrm{^{193m}Hg}$'


# nuclide = r'Hg194'
# LaTeX_RadioNuclide = r'$\mathrm{^{194}Hg}$'


# nuclide = r'Hg195'
# LaTeX_RadioNuclide = r'$\mathrm{^{195}Hg}$'


# nuclide = r'Hg195m'
# LaTeX_RadioNuclide = r'$\mathrm{^{195m}Hg}$'


# nuclide = r'Hg197'
# LaTeX_RadioNuclide = r'$\mathrm{^{197}Hg}$'


# nuclide = r'Hg197m'
# LaTeX_RadioNuclide = r'$\mathrm{^{197m}Hg}$'


# nuclide = r'Hg199m'
# LaTeX_RadioNuclide = r'$\mathrm{^{199m}Hg}$'


# nuclide = r'Hg203'
# LaTeX_RadioNuclide = r'$\mathrm{^{203}Hg}$'


# nuclide = r'Hg205'
# LaTeX_RadioNuclide = r'$\mathrm{^{205}Hg}$'


# nuclide = r'Hg206'
# LaTeX_RadioNuclide = r'$\mathrm{^{206}Hg}$'


# nuclide = r'Hg207'
# LaTeX_RadioNuclide = r'$\mathrm{^{207}Hg}$'


# nuclide = r'Ho150'
# LaTeX_RadioNuclide = r'$\mathrm{^{150}Ho}$'


# nuclide = r'Ho153'
# LaTeX_RadioNuclide = r'$\mathrm{^{153}Ho}$'


# nuclide = r'Ho153m'
# LaTeX_RadioNuclide = r'$\mathrm{^{153m}Ho}$'


# nuclide = r'Ho154'
# LaTeX_RadioNuclide = r'$\mathrm{^{154}Ho}$'


# nuclide = r'Ho154m'
# LaTeX_RadioNuclide = r'$\mathrm{^{154m}Ho}$'


# nuclide = r'Ho155'
# LaTeX_RadioNuclide = r'$\mathrm{^{155}Ho}$'


# nuclide = r'Ho156'
# LaTeX_RadioNuclide = r'$\mathrm{^{156}Ho}$'


# nuclide = r'Ho157'
# LaTeX_RadioNuclide = r'$\mathrm{^{157}Ho}$'


# nuclide = r'Ho159'
# LaTeX_RadioNuclide = r'$\mathrm{^{159}Ho}$'


# nuclide = r'Ho160'
# LaTeX_RadioNuclide = r'$\mathrm{^{160}Ho}$'


# nuclide = r'Ho161'
# LaTeX_RadioNuclide = r'$\mathrm{^{161}Ho}$'


# nuclide = r'Ho162'
# LaTeX_RadioNuclide = r'$\mathrm{^{162}Ho}$'


# nuclide = r'Ho162m'
# LaTeX_RadioNuclide = r'$\mathrm{^{162m}Ho}$'


# nuclide = r'Ho163'
# LaTeX_RadioNuclide = r'$\mathrm{^{163}Ho}$'


# nuclide = r'Ho164'
# LaTeX_RadioNuclide = r'$\mathrm{^{164}Ho}$'


# nuclide = r'Ho164m'
# LaTeX_RadioNuclide = r'$\mathrm{^{164m}Ho}$'


# nuclide = r'Ho166'
# LaTeX_RadioNuclide = r'$\mathrm{^{166}Ho}$'


# nuclide = r'Ho166m'
# LaTeX_RadioNuclide = r'$\mathrm{^{166m}Ho}$'


# nuclide = r'Ho167'
# LaTeX_RadioNuclide = r'$\mathrm{^{167}Ho}$'


# nuclide = r'Ho168'
# LaTeX_RadioNuclide = r'$\mathrm{^{168}Ho}$'


# nuclide = r'Ho168m'
# LaTeX_RadioNuclide = r'$\mathrm{^{168m}Ho}$'


# nuclide = r'Ho170'
# LaTeX_RadioNuclide = r'$\mathrm{^{170}Ho}$'


# nuclide = r'I118'
# LaTeX_RadioNuclide = r'$\mathrm{^{118}I}$'


# nuclide = r'I118m'
# LaTeX_RadioNuclide = r'$\mathrm{^{118m}I}$'


# nuclide = r'I119'
# LaTeX_RadioNuclide = r'$\mathrm{^{119}I}$'


# nuclide = r'I120'
# LaTeX_RadioNuclide = r'$\mathrm{^{120}I}$'


# nuclide = r'I120m'
# LaTeX_RadioNuclide = r'$\mathrm{^{120m}I}$'


# nuclide = r'I121'
# LaTeX_RadioNuclide = r'$\mathrm{^{121}I}$'


# nuclide = r'I122'
# LaTeX_RadioNuclide = r'$\mathrm{^{122}I}$'


# nuclide = r'I123'
# LaTeX_RadioNuclide = r'$\mathrm{^{123}I}$'


# nuclide = r'I124'
# LaTeX_RadioNuclide = r'$\mathrm{^{124}I}$'


# nuclide = r'I125'
# LaTeX_RadioNuclide = r'$\mathrm{^{125}I}$'


# nuclide = r'I126'
# LaTeX_RadioNuclide = r'$\mathrm{^{126}I}$'


# nuclide = r'I128'
# LaTeX_RadioNuclide = r'$\mathrm{^{128}I}$'


# nuclide = r'I129'
# LaTeX_RadioNuclide = r'$\mathrm{^{129}I}$'


# nuclide = r'I130'
# LaTeX_RadioNuclide = r'$\mathrm{^{130}I}$'


# nuclide = r'I130m'
# LaTeX_RadioNuclide = r'$\mathrm{^{130m}I}$'


# nuclide = r'I131'
# LaTeX_RadioNuclide = r'$\mathrm{^{131}I}$'


# nuclide = r'I132'
# LaTeX_RadioNuclide = r'$\mathrm{^{132}I}$'


# nuclide = r'I132m'
# LaTeX_RadioNuclide = r'$\mathrm{^{132m}I}$'


# nuclide = r'I133'
# LaTeX_RadioNuclide = r'$\mathrm{^{133}I}$'


# nuclide = r'I134'
# LaTeX_RadioNuclide = r'$\mathrm{^{134}I}$'


# nuclide = r'I134m'
# LaTeX_RadioNuclide = r'$\mathrm{^{134m}I}$'


# nuclide = r'I135'
# LaTeX_RadioNuclide = r'$\mathrm{^{135}I}$'


# nuclide = r'In103'
# LaTeX_RadioNuclide = r'$\mathrm{^{103}In}$'


# nuclide = r'In105'
# LaTeX_RadioNuclide = r'$\mathrm{^{105}In}$'


# nuclide = r'In106'
# LaTeX_RadioNuclide = r'$\mathrm{^{106}In}$'


# nuclide = r'In106m'
# LaTeX_RadioNuclide = r'$\mathrm{^{106m}In}$'


# nuclide = r'In107'
# LaTeX_RadioNuclide = r'$\mathrm{^{107}In}$'


# nuclide = r'In108'
# LaTeX_RadioNuclide = r'$\mathrm{^{108}In}$'


# nuclide = r'In108m'
# LaTeX_RadioNuclide = r'$\mathrm{^{108m}In}$'


# nuclide = r'In109'
# LaTeX_RadioNuclide = r'$\mathrm{^{109}In}$'


# nuclide = r'In109m'
# LaTeX_RadioNuclide = r'$\mathrm{^{109m}In}$'


# nuclide = r'In110'
# LaTeX_RadioNuclide = r'$\mathrm{^{110}In}$'


# nuclide = r'In110m'
# LaTeX_RadioNuclide = r'$\mathrm{^{110m}In}$'


# nuclide = r'In111'
# LaTeX_RadioNuclide = r'$\mathrm{^{111}In}$'


# nuclide = r'In111m'
# LaTeX_RadioNuclide = r'$\mathrm{^{111m}In}$'


# nuclide = r'In112'
# LaTeX_RadioNuclide = r'$\mathrm{^{112}In}$'


# nuclide = r'In112m'
# LaTeX_RadioNuclide = r'$\mathrm{^{112m}In}$'


# nuclide = r'In113m'
# LaTeX_RadioNuclide = r'$\mathrm{^{113m}In}$'


# nuclide = r'In114'
# LaTeX_RadioNuclide = r'$\mathrm{^{114}In}$'


# nuclide = r'In114m'
# LaTeX_RadioNuclide = r'$\mathrm{^{114m}In}$'


# nuclide = r'In115'
# LaTeX_RadioNuclide = r'$\mathrm{^{115}In}$'


# nuclide = r'In115m'
# LaTeX_RadioNuclide = r'$\mathrm{^{115m}In}$'


# nuclide = r'In116m'
# LaTeX_RadioNuclide = r'$\mathrm{^{116m}In}$'


# nuclide = r'In117'
# LaTeX_RadioNuclide = r'$\mathrm{^{117}In}$'


# nuclide = r'In117m'
# LaTeX_RadioNuclide = r'$\mathrm{^{117m}In}$'


# nuclide = r'In118'
# LaTeX_RadioNuclide = r'$\mathrm{^{118}In}$'


# nuclide = r'In118m'
# LaTeX_RadioNuclide = r'$\mathrm{^{118m}In}$'


# nuclide = r'In119'
# LaTeX_RadioNuclide = r'$\mathrm{^{119}In}$'


# nuclide = r'In119m'
# LaTeX_RadioNuclide = r'$\mathrm{^{119m}In}$'


# nuclide = r'In121'
# LaTeX_RadioNuclide = r'$\mathrm{^{121}In}$'


# nuclide = r'In121m'
# LaTeX_RadioNuclide = r'$\mathrm{^{121m}In}$'


# nuclide = r'Ir180'
# LaTeX_RadioNuclide = r'$\mathrm{^{180}Ir}$'


# nuclide = r'Ir182'
# LaTeX_RadioNuclide = r'$\mathrm{^{182}Ir}$'


# nuclide = r'Ir183'
# LaTeX_RadioNuclide = r'$\mathrm{^{183}Ir}$'


# nuclide = r'Ir184'
# LaTeX_RadioNuclide = r'$\mathrm{^{184}Ir}$'


# nuclide = r'Ir185'
# LaTeX_RadioNuclide = r'$\mathrm{^{185}Ir}$'


# nuclide = r'Ir186'
# LaTeX_RadioNuclide = r'$\mathrm{^{186}Ir}$'


# nuclide = r'Ir186m'
# LaTeX_RadioNuclide = r'$\mathrm{^{186m}Ir}$'


# nuclide = r'Ir187'
# LaTeX_RadioNuclide = r'$\mathrm{^{187}Ir}$'


# nuclide = r'Ir188'
# LaTeX_RadioNuclide = r'$\mathrm{^{188}Ir}$'


# nuclide = r'Ir189'
# LaTeX_RadioNuclide = r'$\mathrm{^{189}Ir}$'


# nuclide = r'Ir190'
# LaTeX_RadioNuclide = r'$\mathrm{^{190}Ir}$'


# nuclide = r'Ir190m'
# LaTeX_RadioNuclide = r'$\mathrm{^{190m}Ir}$'


# nuclide = r'Ir190m2'
# LaTeX_RadioNuclide = r'$\mathrm{^{190m2}Ir}$'


# nuclide = r'Ir191m'
# LaTeX_RadioNuclide = r'$\mathrm{^{191m}Ir}$'


# @@
nuclide = r'Ir192'
LaTeX_RadioNuclide = r'$\mathrm{^{192}Ir}$'


# nuclide = r'Ir192m'
# LaTeX_RadioNuclide = r'$\mathrm{^{192m}Ir}$'


# nuclide = r'Ir192m2'
# LaTeX_RadioNuclide = r'$\mathrm{^{192m2}Ir}$'


# nuclide = r'Ir193m'
# LaTeX_RadioNuclide = r'$\mathrm{^{193m}Ir}$'


# nuclide = r'Ir194'
# LaTeX_RadioNuclide = r'$\mathrm{^{194}Ir}$'


# nuclide = r'Ir194m'
# LaTeX_RadioNuclide = r'$\mathrm{^{194m}Ir}$'


# nuclide = r'Ir195'
# LaTeX_RadioNuclide = r'$\mathrm{^{195}Ir}$'


# nuclide = r'Ir195m'
# LaTeX_RadioNuclide = r'$\mathrm{^{195m}Ir}$'


# nuclide = r'Ir196'
# LaTeX_RadioNuclide = r'$\mathrm{^{196}Ir}$'


# nuclide = r'Ir196m'
# LaTeX_RadioNuclide = r'$\mathrm{^{196m}Ir}$'


# nuclide = r'K38'
# LaTeX_RadioNuclide = r'$\mathrm{^{38}K}$'


# nuclide = r'K40'
# LaTeX_RadioNuclide = r'$\mathrm{_{19}^{40}K}$'


# nuclide = r'K42'
# LaTeX_RadioNuclide = r'$\mathrm{_{19}^{42}K}$'


# nuclide = r'K43'
# LaTeX_RadioNuclide = r'$\mathrm{^{43}K}$'


# nuclide = r'K44'
# LaTeX_RadioNuclide = r'$\mathrm{^{44}K}$'


# nuclide = r'K45'
# LaTeX_RadioNuclide = r'$\mathrm{^{45}K}$'


# nuclide = r'K46'
# LaTeX_RadioNuclide = r'$\mathrm{^{46}K}$'


# nuclide = r'Kr74'
# LaTeX_RadioNuclide = r'$\mathrm{^{74}Kr}$'


# nuclide = r'Kr75'
# LaTeX_RadioNuclide = r'$\mathrm{^{75}Kr}$'


# nuclide = r'Kr76'
# LaTeX_RadioNuclide = r'$\mathrm{^{76}Kr}$'


# nuclide = r'Kr77'
# LaTeX_RadioNuclide = r'$\mathrm{^{77}Kr}$'


# nuclide = r'Kr79'
# LaTeX_RadioNuclide = r'$\mathrm{^{79}Kr}$'


# nuclide = r'Kr81'
# LaTeX_RadioNuclide = r'$\mathrm{^{81}Kr}$'


# nuclide = r'Kr81m'
# LaTeX_RadioNuclide = r'$\mathrm{^{81m}Kr}$'


# nuclide = r'Kr83m'
# LaTeX_RadioNuclide = r'$\mathrm{^{83m}Kr}$'


# nuclide = r'Kr85'
# LaTeX_RadioNuclide = r'$\mathrm{^{85}Kr}$'


# nuclide = r'Kr85m'
# LaTeX_RadioNuclide = r'$\mathrm{^{85m}Kr}$'


# nuclide = r'Kr87'
# LaTeX_RadioNuclide = r'$\mathrm{^{87}Kr}$'


# nuclide = r'Kr88'
# LaTeX_RadioNuclide = r'$\mathrm{^{88}Kr}$'


# nuclide = r'Kr89'
# LaTeX_RadioNuclide = r'$\mathrm{^{89}Kr}$'


# nuclide = r'La128'
# LaTeX_RadioNuclide = r'$\mathrm{^{128}La}$'


# nuclide = r'La129'
# LaTeX_RadioNuclide = r'$\mathrm{^{129}La}$'


# nuclide = r'La130'
# LaTeX_RadioNuclide = r'$\mathrm{^{130}La}$'


# nuclide = r'La131'
# LaTeX_RadioNuclide = r'$\mathrm{^{131}La}$'


# nuclide = r'La132'
# LaTeX_RadioNuclide = r'$\mathrm{^{132}La}$'


# nuclide = r'La132m'
# LaTeX_RadioNuclide = r'$\mathrm{^{132m}La}$'


# nuclide = r'La133'
# LaTeX_RadioNuclide = r'$\mathrm{^{133}La}$'


# nuclide = r'La134'
# LaTeX_RadioNuclide = r'$\mathrm{^{134}La}$'


# nuclide = r'La135'
# LaTeX_RadioNuclide = r'$\mathrm{^{135}La}$'


# nuclide = r'La136'
# LaTeX_RadioNuclide = r'$\mathrm{^{136}La}$'


# nuclide = r'La137'
# LaTeX_RadioNuclide = r'$\mathrm{^{137}La}$'


# nuclide = r'La138'
# LaTeX_RadioNuclide = r'$\mathrm{^{138}La}$'


# nuclide = r'La140'
# LaTeX_RadioNuclide = r'$\mathrm{^{140}La}$'


# nuclide = r'La141'
# LaTeX_RadioNuclide = r'$\mathrm{^{141}La}$'


# nuclide = r'La142'
# LaTeX_RadioNuclide = r'$\mathrm{^{142}La}$'


# nuclide = r'La143'
# LaTeX_RadioNuclide = r'$\mathrm{^{143}La}$'


# nuclide = r'Lu165'
# LaTeX_RadioNuclide = r'$\mathrm{^{165}Lu}$'


# nuclide = r'Lu167'
# LaTeX_RadioNuclide = r'$\mathrm{^{167}Lu}$'


# nuclide = r'Lu169'
# LaTeX_RadioNuclide = r'$\mathrm{^{169}Lu}$'


# nuclide = r'Lu169m'
# LaTeX_RadioNuclide = r'$\mathrm{^{169m}Lu}$'


# nuclide = r'Lu170'
# LaTeX_RadioNuclide = r'$\mathrm{^{170}Lu}$'


# nuclide = r'Lu171'
# LaTeX_RadioNuclide = r'$\mathrm{^{171}Lu}$'


# nuclide = r'Lu171m'
# LaTeX_RadioNuclide = r'$\mathrm{^{171m}Lu}$'


# nuclide = r'Lu172'
# LaTeX_RadioNuclide = r'$\mathrm{^{172}Lu}$'


# nuclide = r'Lu172m'
# LaTeX_RadioNuclide = r'$\mathrm{^{172m}Lu}$'


# nuclide = r'Lu173'
# LaTeX_RadioNuclide = r'$\mathrm{^{173}Lu}$'


# nuclide = r'Lu174'
# LaTeX_RadioNuclide = r'$\mathrm{^{174}Lu}$'


# nuclide = r'Lu174m'
# LaTeX_RadioNuclide = r'$\mathrm{^{174m}Lu}$'


# nuclide = r'Lu176'
# LaTeX_RadioNuclide = r'$\mathrm{^{176}Lu}$'


# nuclide = r'Lu176m'
# LaTeX_RadioNuclide = r'$\mathrm{^{176m}Lu}$'


# nuclide = r'Lu177'
# LaTeX_RadioNuclide = r'$\mathrm{^{177}Lu}$'


# nuclide = r'Lu177m'
# LaTeX_RadioNuclide = r'$\mathrm{^{177m}Lu}$'


# nuclide = r'Lu178'
# LaTeX_RadioNuclide = r'$\mathrm{^{178}Lu}$'


# nuclide = r'Lu178m'
# LaTeX_RadioNuclide = r'$\mathrm{^{178m}Lu}$'


# nuclide = r'Lu179'
# LaTeX_RadioNuclide = r'$\mathrm{^{179}Lu}$'


# nuclide = r'Lu180'
# LaTeX_RadioNuclide = r'$\mathrm{^{180}Lu}$'


# nuclide = r'Lu181'
# LaTeX_RadioNuclide = r'$\mathrm{^{181}Lu}$'


# nuclide = r'Mg27'
# LaTeX_RadioNuclide = r'$\mathrm{^{27}Mg}$'


# nuclide = r'Mg28'
# LaTeX_RadioNuclide = r'$\mathrm{^{28}Mg}$'


# nuclide = r'Mn50m'
# LaTeX_RadioNuclide = r'$\mathrm{^{50m}Mn}$'


# nuclide = r'Mn51'
# LaTeX_RadioNuclide = r'$\mathrm{^{51}Mn}$'


# nuclide = r'Mn52'
# LaTeX_RadioNuclide = r'$\mathrm{^{52}Mn}$'


# nuclide = r'Mn52m'
# LaTeX_RadioNuclide = r'$\mathrm{^{52m}Mn}$'


# nuclide = r'Mn53'
# LaTeX_RadioNuclide = r'$\mathrm{_{25}^{53}Mn}$'


# nuclide = r'Mn54'
# LaTeX_RadioNuclide = r'$\mathrm{_{25}^{54}Mn}$'


# nuclide = r'Mn56'
# LaTeX_RadioNuclide = r'$\mathrm{_{25}^{56}Mn}$'


# nuclide = r'Mn57'
# LaTeX_RadioNuclide = r'$\mathrm{^{57}Mn}$'


# nuclide = r'Mn58m'
# LaTeX_RadioNuclide = r'$\mathrm{^{58m}Mn}$'


# nuclide = r'Mo89'
# LaTeX_RadioNuclide = r'$\mathrm{^{89}Mo}$'


# nuclide = r'Mo90'
# LaTeX_RadioNuclide = r'$\mathrm{^{90}Mo}$'


# nuclide = r'Mo91'
# LaTeX_RadioNuclide = r'$\mathrm{^{91}Mo}$'


# nuclide = r'Mo91m'
# LaTeX_RadioNuclide = r'$\mathrm{^{91m}Mo}$'


# nuclide = r'Mo93'
# LaTeX_RadioNuclide = r'$\mathrm{_{42}^{93}Mo}$'


# nuclide = r'Mo93m'
# LaTeX_RadioNuclide = r'$\mathrm{^{93m}Mo}$'


# nuclide = r'Mo99'
# LaTeX_RadioNuclide = r'$\mathrm{_{42}^{99}Mo}$'


# nuclide = r'Mo99_Tc99m'
# LaTeX_RadioNuclide = r'$\mathrm{^{99}Mo} + \mathrm{^{99m}Tc}$'


# nuclide = r'Mo101'
# LaTeX_RadioNuclide = r'$\mathrm{^{101}Mo}$'


# nuclide = r'Mo102'
# LaTeX_RadioNuclide = r'$\mathrm{^{102}Mo}$'


# nuclide = r'N13'
# LaTeX_RadioNuclide = r'$\mathrm{^{13}N}$'


# nuclide = r'N16'
# LaTeX_RadioNuclide = r'$\mathrm{^{16}N}$'


# nuclide = r'Na22'
# LaTeX_RadioNuclide = r'$\mathrm{_{11}^{22}Na}$'


# nuclide = r'Na24'
# LaTeX_RadioNuclide = r'$\mathrm{_{11}^{24}Na}$'


# nuclide = r'Nb87'
# LaTeX_RadioNuclide = r'$\mathrm{^{87}Nb}$'


# nuclide = r'Nb88'
# LaTeX_RadioNuclide = r'$\mathrm{^{88}Nb}$'


# nuclide = r'Nb88m'
# LaTeX_RadioNuclide = r'$\mathrm{^{88m}Nb}$'


# nuclide = r'Nb89'
# LaTeX_RadioNuclide = r'$\mathrm{^{89}Nb}$'


# nuclide = r'Nb89m'
# LaTeX_RadioNuclide = r'$\mathrm{^{89m}Nb}$'


# nuclide = r'Nb90'
# LaTeX_RadioNuclide = r'$\mathrm{^{90}Nb}$'


# nuclide = r'Nb91'
# LaTeX_RadioNuclide = r'$\mathrm{^{91}Nb}$'


# nuclide = r'Nb91m'
# LaTeX_RadioNuclide = r'$\mathrm{^{91m}Nb}$'


# nuclide = r'Nb92'
# LaTeX_RadioNuclide = r'$\mathrm{^{92}Nb}$'


# nuclide = r'Nb92m'
# LaTeX_RadioNuclide = r'$\mathrm{^{92m}Nb}$'


# nuclide = r'Nb93m'
# LaTeX_RadioNuclide = r'$\mathrm{^{93m}Nb}$'


# nuclide = r'Nb94'
# LaTeX_RadioNuclide = r'$\mathrm{_{41}^{94}Nb}$'


# nuclide = r'Nb94m'
# LaTeX_RadioNuclide = r'$\mathrm{^{94m}Nb}$'


# nuclide = r'Nb95'
# LaTeX_RadioNuclide = r'$\mathrm{_{41}^{95}Nb}$'


# nuclide = r'Nb95m'
# LaTeX_RadioNuclide = r'$\mathrm{^{95m}Nb}$'


# nuclide = r'Nb96'
# LaTeX_RadioNuclide = r'$\mathrm{^{96}Nb}$'


# nuclide = r'Nb97'
# LaTeX_RadioNuclide = r'$\mathrm{^{97}Nb}$'


# nuclide = r'Nb98m'
# LaTeX_RadioNuclide = r'$\mathrm{^{98m}Nb}$'


# nuclide = r'Nb99'
# LaTeX_RadioNuclide = r'$\mathrm{^{99}Nb}$'


# nuclide = r'Nb99m'
# LaTeX_RadioNuclide = r'$\mathrm{^{99m}Nb}$'


# nuclide = r'Nd134'
# LaTeX_RadioNuclide = r'$\mathrm{^{134}Nd}$'


# nuclide = r'Nd135'
# LaTeX_RadioNuclide = r'$\mathrm{^{135}Nd}$'


# nuclide = r'Nd136'
# LaTeX_RadioNuclide = r'$\mathrm{^{136}Nd}$'


# nuclide = r'Nd137'
# LaTeX_RadioNuclide = r'$\mathrm{^{137}Nd}$'


# nuclide = r'Nd138'
# LaTeX_RadioNuclide = r'$\mathrm{^{138}Nd}$'


# nuclide = r'Nd139'
# LaTeX_RadioNuclide = r'$\mathrm{^{139}Nd}$'


# nuclide = r'Nd139m'
# LaTeX_RadioNuclide = r'$\mathrm{^{139m}Nd}$'


# nuclide = r'Nd140'
# LaTeX_RadioNuclide = r'$\mathrm{^{140}Nd}$'


# nuclide = r'Nd141'
# LaTeX_RadioNuclide = r'$\mathrm{^{141}Nd}$'


# nuclide = r'Nd141m'
# LaTeX_RadioNuclide = r'$\mathrm{^{141m}Nd}$'


# nuclide = r'Nd144'
# LaTeX_RadioNuclide = r'$\mathrm{^{144}Nd}$'


# nuclide = r'Nd147'
# LaTeX_RadioNuclide = r'$\mathrm{^{147}Nd}$'


# nuclide = r'Nd149'
# LaTeX_RadioNuclide = r'$\mathrm{^{149}Nd}$'


# nuclide = r'Nd151'
# LaTeX_RadioNuclide = r'$\mathrm{^{151}Nd}$'


# nuclide = r'Nd152'
# LaTeX_RadioNuclide = r'$\mathrm{^{152}Nd}$'


# nuclide = r'Ne19'
# LaTeX_RadioNuclide = r'$\mathrm{^{19}Ne}$'


# nuclide = r'Ne24'
# LaTeX_RadioNuclide = r'$\mathrm{^{24}Ne}$'


# nuclide = r'Ni56'
# LaTeX_RadioNuclide = r'$\mathrm{^{56}Ni}$'


# nuclide = r'Ni57'
# LaTeX_RadioNuclide = r'$\mathrm{^{57}Ni}$'


# nuclide = r'Ni59'
# LaTeX_RadioNuclide = r'$\mathrm{_{28}^{59}Ni}$'


# nuclide = r'Ni63'
# LaTeX_RadioNuclide = r'$\mathrm{_{28}^{63}Ni}$'


# nuclide = r'Ni65'
# LaTeX_RadioNuclide = r'$\mathrm{^{65}Ni}$'


# nuclide = r'Ni66'
# LaTeX_RadioNuclide = r'$\mathrm{^{66}Ni}$'


# nuclide = r'Np232'
# LaTeX_RadioNuclide = r'$\mathrm{^{232}Np}$'


# nuclide = r'Np233'
# LaTeX_RadioNuclide = r'$\mathrm{^{233}Np}$'


# nuclide = r'Np234'
# LaTeX_RadioNuclide = r'$\mathrm{^{234}Np}$'


# nuclide = r'Np235'
# LaTeX_RadioNuclide = r'$\mathrm{^{235}Np}$'


# nuclide = r'Np236'
# LaTeX_RadioNuclide = r'$\mathrm{^{236}Np}$'


# nuclide = r'Np236m'
# LaTeX_RadioNuclide = r'$\mathrm{^{236m}Np}$'


# nuclide = r'Np237'
# LaTeX_RadioNuclide = r'$\mathrm{^{237}Np}$'


# nuclide = r'Np238'
# LaTeX_RadioNuclide = r'$\mathrm{^{238}Np}$'


# nuclide = r'Np239'
# LaTeX_RadioNuclide = r'$\mathrm{^{239}Np}$'


# nuclide = r'Np240'
# LaTeX_RadioNuclide = r'$\mathrm{^{240}Np}$'


# nuclide = r'Np240m'
# LaTeX_RadioNuclide = r'$\mathrm{^{240m}Np}$'


# nuclide = r'Np241'
# LaTeX_RadioNuclide = r'$\mathrm{^{241}Np}$'


# nuclide = r'Np242'
# LaTeX_RadioNuclide = r'$\mathrm{^{242}Np}$'


# nuclide = r'Np242m'
# LaTeX_RadioNuclide = r'$\mathrm{^{242m}Np}$'


# nuclide = r'O14'
# LaTeX_RadioNuclide = r'$\mathrm{^{14}O}$'


# nuclide = r'O15'
# LaTeX_RadioNuclide = r'$\mathrm{^{15}O}$'


# nuclide = r'O19'
# LaTeX_RadioNuclide = r'$\mathrm{^{19}O}$'


# nuclide = r'Os180'
# LaTeX_RadioNuclide = r'$\mathrm{^{180}Os}$'


# nuclide = r'Os181'
# LaTeX_RadioNuclide = r'$\mathrm{^{181}Os}$'


# nuclide = r'Os182'
# LaTeX_RadioNuclide = r'$\mathrm{^{182}Os}$'


# nuclide = r'Os183'
# LaTeX_RadioNuclide = r'$\mathrm{^{183}Os}$'


# nuclide = r'Os183m'
# LaTeX_RadioNuclide = r'$\mathrm{^{183m}Os}$'


# nuclide = r'Os185'
# LaTeX_RadioNuclide = r'$\mathrm{^{185}Os}$'


# nuclide = r'Os186'
# LaTeX_RadioNuclide = r'$\mathrm{^{186}Os}$'


# nuclide = r'Os189m'
# LaTeX_RadioNuclide = r'$\mathrm{^{189m}Os}$'


# nuclide = r'Os190m'
# LaTeX_RadioNuclide = r'$\mathrm{^{190m}Os}$'


# nuclide = r'Os191'
# LaTeX_RadioNuclide = r'$\mathrm{^{191}Os}$'


# nuclide = r'Os191m'
# LaTeX_RadioNuclide = r'$\mathrm{^{191m}Os}$'


# nuclide = r'Os193'
# LaTeX_RadioNuclide = r'$\mathrm{^{193}Os}$'


# nuclide = r'Os194'
# LaTeX_RadioNuclide = r'$\mathrm{^{194}Os}$'


# nuclide = r'Os196'
# LaTeX_RadioNuclide = r'$\mathrm{^{196}Os}$'


# nuclide = r'P30'
# LaTeX_RadioNuclide = r'$\mathrm{^{30}P}$'


# nuclide = r'P32'
# LaTeX_RadioNuclide = r'$\mathrm{_{15}^{32}P}$'


# nuclide = r'P33'
# LaTeX_RadioNuclide = r'$\mathrm{_{15}^{33}P}$'


# nuclide = r'Pa227'
# LaTeX_RadioNuclide = r'$\mathrm{^{227}Pa}$'


# nuclide = r'Pa228'
# LaTeX_RadioNuclide = r'$\mathrm{^{228}Pa}$'


# nuclide = r'Pa229'
# LaTeX_RadioNuclide = r'$\mathrm{^{229}Pa}$'


# nuclide = r'Pa230'
# LaTeX_RadioNuclide = r'$\mathrm{^{230}Pa}$'


# nuclide = r'Pa231'
# LaTeX_RadioNuclide = r'$\mathrm{^{231}Pa}$'


# nuclide = r'Pa232'
# LaTeX_RadioNuclide = r'$\mathrm{^{232}Pa}$'


# nuclide = r'Pa233'
# LaTeX_RadioNuclide = r'$\mathrm{^{233}Pa}$'


# nuclide = r'Pa234'
# LaTeX_RadioNuclide = r'$\mathrm{^{234}Pa}$'


# nuclide = r'Pa234m'
# LaTeX_RadioNuclide = r'$\mathrm{^{234m}Pa}$'


# nuclide = r'Pa235'
# LaTeX_RadioNuclide = r'$\mathrm{^{235}Pa}$'


# nuclide = r'Pa236'
# LaTeX_RadioNuclide = r'$\mathrm{^{236}Pa}$'


# nuclide = r'Pa237'
# LaTeX_RadioNuclide = r'$\mathrm{^{237}Pa}$'


# nuclide = r'Pb194'
# LaTeX_RadioNuclide = r'$\mathrm{^{194}Pb}$'


# nuclide = r'Pb195m'
# LaTeX_RadioNuclide = r'$\mathrm{^{195m}Pb}$'


# nuclide = r'Pb196'
# LaTeX_RadioNuclide = r'$\mathrm{^{196}Pb}$'


# nuclide = r'Pb197'
# LaTeX_RadioNuclide = r'$\mathrm{^{197}Pb}$'


# nuclide = r'Pb197m'
# LaTeX_RadioNuclide = r'$\mathrm{^{197m}Pb}$'


# nuclide = r'Pb198'
# LaTeX_RadioNuclide = r'$\mathrm{^{198}Pb}$'


# nuclide = r'Pb199'
# LaTeX_RadioNuclide = r'$\mathrm{^{199}Pb}$'


# nuclide = r'Pb200'
# LaTeX_RadioNuclide = r'$\mathrm{^{200}Pb}$'


# nuclide = r'Pb201'
# LaTeX_RadioNuclide = r'$\mathrm{^{201}Pb}$'


# nuclide = r'Pb201m'
# LaTeX_RadioNuclide = r'$\mathrm{^{201m}Pb}$'


# nuclide = r'Pb202'
# LaTeX_RadioNuclide = r'$\mathrm{^{202}Pb}$'


# nuclide = r'Pb202m'
# LaTeX_RadioNuclide = r'$\mathrm{^{202m}Pb}$'


# nuclide = r'Pb203'
# LaTeX_RadioNuclide = r'$\mathrm{^{203}Pb}$'


# nuclide = r'Pb204m'
# LaTeX_RadioNuclide = r'$\mathrm{^{204m}Pb}$'


# nuclide = r'Pb205'
# LaTeX_RadioNuclide = r'$\mathrm{^{205}Pb}$'


# nuclide = r'Pb209'
# LaTeX_RadioNuclide = r'$\mathrm{^{209}Pb}$'


# nuclide = r'Pb210'
# LaTeX_RadioNuclide = r'$\mathrm{^{210}Pb}$'


# nuclide = r'Pb211'
# LaTeX_RadioNuclide = r'$\mathrm{^{211}Pb}$'


# nuclide = r'Pb212'
# LaTeX_RadioNuclide = r'$\mathrm{^{212}Pb}$'


# nuclide = r'Pb214'
# LaTeX_RadioNuclide = r'$\mathrm{^{214}Pb}$'


# nuclide = r'Pd96'
# LaTeX_RadioNuclide = r'$\mathrm{^{96}Pd}$'


# nuclide = r'Pd97'
# LaTeX_RadioNuclide = r'$\mathrm{^{97}Pd}$'


# nuclide = r'Pd98'
# LaTeX_RadioNuclide = r'$\mathrm{^{98}Pd}$'


# nuclide = r'Pd99'
# LaTeX_RadioNuclide = r'$\mathrm{^{99}Pd}$'


# nuclide = r'Pd100'
# LaTeX_RadioNuclide = r'$\mathrm{^{100}Pd}$'


# nuclide = r'Pd101'
# LaTeX_RadioNuclide = r'$\mathrm{^{101}Pd}$'


# nuclide = r'Pd103'
# LaTeX_RadioNuclide = r'$\mathrm{^{103}Pd}$'


# nuclide = r'Pd103_Rh103m'
# LaTeX_RadioNuclide = r'$\mathrm{^{103}Pd} + \mathrm{^{103m}Rh}$'


# nuclide = r'Pd107'
# LaTeX_RadioNuclide = r'$\mathrm{^{107}Pd}$'


# nuclide = r'Pd109'
# LaTeX_RadioNuclide = r'$\mathrm{^{109}Pd}$'


# nuclide = r'Pd109m'
# LaTeX_RadioNuclide = r'$\mathrm{^{109m}Pd}$'


# nuclide = r'Pd111'
# LaTeX_RadioNuclide = r'$\mathrm{^{111}Pd}$'


# nuclide = r'Pd112'
# LaTeX_RadioNuclide = r'$\mathrm{^{112}Pd}$'


# nuclide = r'Pd114'
# LaTeX_RadioNuclide = r'$\mathrm{^{114}Pd}$'


# nuclide = r'Pm136'
# LaTeX_RadioNuclide = r'$\mathrm{^{136}Pm}$'


# nuclide = r'Pm137m'
# LaTeX_RadioNuclide = r'$\mathrm{^{137m}Pm}$'


# nuclide = r'Pm139'
# LaTeX_RadioNuclide = r'$\mathrm{^{139}Pm}$'


# nuclide = r'Pm140'
# LaTeX_RadioNuclide = r'$\mathrm{^{140}Pm}$'


# nuclide = r'Pm140m'
# LaTeX_RadioNuclide = r'$\mathrm{^{140m}Pm}$'


# nuclide = r'Pm141'
# LaTeX_RadioNuclide = r'$\mathrm{^{141}Pm}$'


# nuclide = r'Pm142'
# LaTeX_RadioNuclide = r'$\mathrm{^{142}Pm}$'


# nuclide = r'Pm143'
# LaTeX_RadioNuclide = r'$\mathrm{^{143}Pm}$'


# nuclide = r'Pm144'
# LaTeX_RadioNuclide = r'$\mathrm{^{144}Pm}$'


# nuclide = r'Pm145'
# LaTeX_RadioNuclide = r'$\mathrm{^{145}Pm}$'


# nuclide = r'Pm146'
# LaTeX_RadioNuclide = r'$\mathrm{^{146}Pm}$'


# nuclide = r'Pm147'
# LaTeX_RadioNuclide = r'$\mathrm{^{147}Pm}$'


# nuclide = r'Pm148'
# LaTeX_RadioNuclide = r'$\mathrm{^{148}Pm}$'


# nuclide = r'Pm148m'
# LaTeX_RadioNuclide = r'$\mathrm{^{148m}Pm}$'


# nuclide = r'Pm149'
# LaTeX_RadioNuclide = r'$\mathrm{^{149}Pm}$'


# nuclide = r'Pm150'
# LaTeX_RadioNuclide = r'$\mathrm{^{150}Pm}$'


# nuclide = r'Pm151'
# LaTeX_RadioNuclide = r'$\mathrm{^{151}Pm}$'


# nuclide = r'Pm152'
# LaTeX_RadioNuclide = r'$\mathrm{^{152}Pm}$'


# nuclide = r'Pm152m'
# LaTeX_RadioNuclide = r'$\mathrm{^{152m}Pm}$'


# nuclide = r'Pm153'
# LaTeX_RadioNuclide = r'$\mathrm{^{153}Pm}$'


# nuclide = r'Pm154'
# LaTeX_RadioNuclide = r'$\mathrm{^{154}Pm}$'


# nuclide = r'Pm154m'
# LaTeX_RadioNuclide = r'$\mathrm{^{154m}Pm}$'


# nuclide = r'Po203'
# LaTeX_RadioNuclide = r'$\mathrm{^{203}Po}$'


# nuclide = r'Po204'
# LaTeX_RadioNuclide = r'$\mathrm{^{204}Po}$'


# nuclide = r'Po205'
# LaTeX_RadioNuclide = r'$\mathrm{^{205}Po}$'


# nuclide = r'Po206'
# LaTeX_RadioNuclide = r'$\mathrm{^{206}Po}$'


# nuclide = r'Po207'
# LaTeX_RadioNuclide = r'$\mathrm{^{207}Po}$'


# nuclide = r'Po208'
# LaTeX_RadioNuclide = r'$\mathrm{^{208}Po}$'


# nuclide = r'Po209'
# LaTeX_RadioNuclide = r'$\mathrm{^{209}Po}$'


# nuclide = r'Po210'
# LaTeX_RadioNuclide = r'$\mathrm{^{210}Po}$'


# nuclide = r'Po211'
# LaTeX_RadioNuclide = r'$\mathrm{^{211}Po}$'


# nuclide = r'Po212'
# LaTeX_RadioNuclide = r'$\mathrm{^{212}Po}$'


# nuclide = r'Po212m'
# LaTeX_RadioNuclide = r'$\mathrm{^{212m}Po}$'


# nuclide = r'Po213'
# LaTeX_RadioNuclide = r'$\mathrm{^{213}Po}$'


# nuclide = r'Po214'
# LaTeX_RadioNuclide = r'$\mathrm{^{214}Po}$'


# nuclide = r'Po215'
# LaTeX_RadioNuclide = r'$\mathrm{^{215}Po}$'


# nuclide = r'Po216'
# LaTeX_RadioNuclide = r'$\mathrm{^{216}Po}$'


# nuclide = r'Po218'
# LaTeX_RadioNuclide = r'$\mathrm{^{218}Po}$'


# nuclide = r'Pr134'
# LaTeX_RadioNuclide = r'$\mathrm{^{134}Pr}$'


# nuclide = r'Pr134m'
# LaTeX_RadioNuclide = r'$\mathrm{^{134m}Pr}$'


# nuclide = r'Pr135'
# LaTeX_RadioNuclide = r'$\mathrm{^{135}Pr}$'


# nuclide = r'Pr136'
# LaTeX_RadioNuclide = r'$\mathrm{^{136}Pr}$'


# nuclide = r'Pr137'
# LaTeX_RadioNuclide = r'$\mathrm{^{137}Pr}$'


# nuclide = r'Pr138'
# LaTeX_RadioNuclide = r'$\mathrm{^{138}Pr}$'


# nuclide = r'Pr138m'
# LaTeX_RadioNuclide = r'$\mathrm{^{138m}Pr}$'


# nuclide = r'Pr139'
# LaTeX_RadioNuclide = r'$\mathrm{^{139}Pr}$'


# nuclide = r'Pr140'
# LaTeX_RadioNuclide = r'$\mathrm{^{140}Pr}$'


# nuclide = r'Pr142'
# LaTeX_RadioNuclide = r'$\mathrm{^{142}Pr}$'


# nuclide = r'Pr142m'
# LaTeX_RadioNuclide = r'$\mathrm{^{142m}Pr}$'


# nuclide = r'Pr143'
# LaTeX_RadioNuclide = r'$\mathrm{^{143}Pr}$'


# nuclide = r'Pr144'
# LaTeX_RadioNuclide = r'$\mathrm{^{144}Pr}$'


# nuclide = r'Pr144m'
# LaTeX_RadioNuclide = r'$\mathrm{^{144m}Pr}$'


# nuclide = r'Pr145'
# LaTeX_RadioNuclide = r'$\mathrm{^{145}Pr}$'


# nuclide = r'Pr146'
# LaTeX_RadioNuclide = r'$\mathrm{^{146}Pr}$'


# nuclide = r'Pr147'
# LaTeX_RadioNuclide = r'$\mathrm{^{147}Pr}$'


# nuclide = r'Pr148'
# LaTeX_RadioNuclide = r'$\mathrm{^{148}Pr}$'


# nuclide = r'Pr148m'
# LaTeX_RadioNuclide = r'$\mathrm{^{148m}Pr}$'


# nuclide = r'Pt184'
# LaTeX_RadioNuclide = r'$\mathrm{^{184}Pt}$'


# nuclide = r'Pt186'
# LaTeX_RadioNuclide = r'$\mathrm{^{186}Pt}$'


# nuclide = r'Pt187'
# LaTeX_RadioNuclide = r'$\mathrm{^{187}Pt}$'


# nuclide = r'Pt188'
# LaTeX_RadioNuclide = r'$\mathrm{^{188}Pt}$'


# nuclide = r'Pt189'
# LaTeX_RadioNuclide = r'$\mathrm{^{189}Pt}$'


# nuclide = r'Pt190'
# LaTeX_RadioNuclide = r'$\mathrm{^{190}Pt}$'


# nuclide = r'Pt191'
# LaTeX_RadioNuclide = r'$\mathrm{^{191}Pt}$'


# nuclide = r'Pt193'
# LaTeX_RadioNuclide = r'$\mathrm{^{193}Pt}$'


# nuclide = r'Pt193m'
# LaTeX_RadioNuclide = r'$\mathrm{^{193m}Pt}$'


# nuclide = r'Pt195m'
# LaTeX_RadioNuclide = r'$\mathrm{^{195m}Pt}$'


# nuclide = r'Pt197'
# LaTeX_RadioNuclide = r'$\mathrm{^{197}Pt}$'


# nuclide = r'Pt197m'
# LaTeX_RadioNuclide = r'$\mathrm{^{197m}Pt}$'


# nuclide = r'Pt199'
# LaTeX_RadioNuclide = r'$\mathrm{^{199}Pt}$'


# nuclide = r'Pt200'
# LaTeX_RadioNuclide = r'$\mathrm{^{200}Pt}$'


# nuclide = r'Pt202'
# LaTeX_RadioNuclide = r'$\mathrm{^{202}Pt}$'


# nuclide = r'Pu232'
# LaTeX_RadioNuclide = r'$\mathrm{^{232}Pu}$'


# nuclide = r'Pu234'
# LaTeX_RadioNuclide = r'$\mathrm{^{234}Pu}$'


# nuclide = r'Pu235'
# LaTeX_RadioNuclide = r'$\mathrm{^{235}Pu}$'


# nuclide = r'Pu236'
# LaTeX_RadioNuclide = r'$\mathrm{^{236}Pu}$'


# nuclide = r'Pu237'
# LaTeX_RadioNuclide = r'$\mathrm{^{237}Pu}$'


# nuclide = r'Pu238'
# LaTeX_RadioNuclide = r'$\mathrm{^{238}Pu}$'


# nuclide = r'Pu239'
# LaTeX_RadioNuclide = r'$\mathrm{^{239}Pu}$'


# nuclide = r'Pu240'
# LaTeX_RadioNuclide = r'$\mathrm{^{240}Pu}$'


# nuclide = r'Pu241'
# LaTeX_RadioNuclide = r'$\mathrm{^{241}Pu}$'


# nuclide = r'Pu242'
# LaTeX_RadioNuclide = r'$\mathrm{^{242}Pu}$'


# nuclide = r'Pu243'
# LaTeX_RadioNuclide = r'$\mathrm{^{243}Pu}$'


# nuclide = r'Pu244'
# LaTeX_RadioNuclide = r'$\mathrm{^{244}Pu}$'


# nuclide = r'Pu245'
# LaTeX_RadioNuclide = r'$\mathrm{^{245}Pu}$'


# nuclide = r'Pu246'
# LaTeX_RadioNuclide = r'$\mathrm{^{246}Pu}$'


# nuclide = r'Ra219'
# LaTeX_RadioNuclide = r'$\mathrm{^{219}Ra}$'


# nuclide = r'Ra220'
# LaTeX_RadioNuclide = r'$\mathrm{^{220}Ra}$'


# nuclide = r'Ra221'
# LaTeX_RadioNuclide = r'$\mathrm{^{221}Ra}$'


# nuclide = r'Ra222'
# LaTeX_RadioNuclide = r'$\mathrm{^{222}Ra}$'


# nuclide = r'Ra223'
# LaTeX_RadioNuclide = r'$\mathrm{^{223}Ra}$'


# nuclide = r'Ra224'
# LaTeX_RadioNuclide = r'$\mathrm{^{224}Ra}$'


# nuclide = r'Ra225'
# LaTeX_RadioNuclide = r'$\mathrm{^{225}Ra}$'


# nuclide = r'Ra226'
# LaTeX_RadioNuclide = r'$\mathrm{^{226}Ra}$'


# nuclide = r'Ra227'
# LaTeX_RadioNuclide = r'$\mathrm{^{227}Ra}$'


# nuclide = r'Ra228'
# LaTeX_RadioNuclide = r'$\mathrm{^{228}Ra}$'


# nuclide = r'Ra230'
# LaTeX_RadioNuclide = r'$\mathrm{^{230}Ra}$'


# nuclide = r'Rb77'
# LaTeX_RadioNuclide = r'$\mathrm{^{77}Rb}$'


# nuclide = r'Rb78'
# LaTeX_RadioNuclide = r'$\mathrm{^{78}Rb}$'


# nuclide = r'Rb78m'
# LaTeX_RadioNuclide = r'$\mathrm{^{78m}Rb}$'


# nuclide = r'Rb79'
# LaTeX_RadioNuclide = r'$\mathrm{^{79}Rb}$'


# nuclide = r'Rb80'
# LaTeX_RadioNuclide = r'$\mathrm{^{80}Rb}$'


# nuclide = r'Rb81'
# LaTeX_RadioNuclide = r'$\mathrm{^{81}Rb}$'


# nuclide = r'Rb81m'
# LaTeX_RadioNuclide = r'$\mathrm{^{81m}Rb}$'


# nuclide = r'Rb82'
# LaTeX_RadioNuclide = r'$\mathrm{_{37}^{82}Rb}$'


# nuclide = r'Rb82m'
# LaTeX_RadioNuclide = r'$\mathrm{^{82m}Rb}$'


# nuclide = r'Rb83'
# LaTeX_RadioNuclide = r'$\mathrm{^{83}Rb}$'


# nuclide = r'Rb84'
# LaTeX_RadioNuclide = r'$\mathrm{^{84}Rb}$'


# nuclide = r'Rb84m'
# LaTeX_RadioNuclide = r'$\mathrm{^{84m}Rb}$'


# nuclide = r'Rb86'
# LaTeX_RadioNuclide = r'$\mathrm{^{86}Rb}$'


# nuclide = r'Rb86m'
# LaTeX_RadioNuclide = r'$\mathrm{^{86m}Rb}$'


# nuclide = r'Rb87'
# LaTeX_RadioNuclide = r'$\mathrm{^{87}Rb}$'


# nuclide = r'Rb88'
# LaTeX_RadioNuclide = r'$\mathrm{^{88}Rb}$'


# nuclide = r'Rb89'
# LaTeX_RadioNuclide = r'$\mathrm{^{89}Rb}$'


# nuclide = r'Rb90'
# LaTeX_RadioNuclide = r'$\mathrm{^{90}Rb}$'


# nuclide = r'Rb90m'
# LaTeX_RadioNuclide = r'$\mathrm{^{90m}Rb}$'


# nuclide = r'Re178'
# LaTeX_RadioNuclide = r'$\mathrm{^{178}Re}$'


# nuclide = r'Re179'
# LaTeX_RadioNuclide = r'$\mathrm{^{179}Re}$'


# nuclide = r'Re180'
# LaTeX_RadioNuclide = r'$\mathrm{^{180}Re}$'


# nuclide = r'Re181'
# LaTeX_RadioNuclide = r'$\mathrm{^{181}Re}$'


# nuclide = r'Re182'
# LaTeX_RadioNuclide = r'$\mathrm{^{182}Re}$'


# nuclide = r'Re182m'
# LaTeX_RadioNuclide = r'$\mathrm{^{182m}Re}$'


# nuclide = r'Re183'
# LaTeX_RadioNuclide = r'$\mathrm{^{183}Re}$'


# nuclide = r'Re184'
# LaTeX_RadioNuclide = r'$\mathrm{^{184}Re}$'


# nuclide = r'Re184m'
# LaTeX_RadioNuclide = r'$\mathrm{^{184m}Re}$'


# nuclide = r'Re186'
# LaTeX_RadioNuclide = r'$\mathrm{^{186}Re}$'


# nuclide = r'Re186m'
# LaTeX_RadioNuclide = r'$\mathrm{^{186m}Re}$'


# nuclide = r'Re187'
# LaTeX_RadioNuclide = r'$\mathrm{^{187}Re}$'


# nuclide = r'Re188'
# LaTeX_RadioNuclide = r'$\mathrm{^{188}Re}$'


# nuclide = r'Re188m'
# LaTeX_RadioNuclide = r'$\mathrm{^{188m}Re}$'


# nuclide = r'Re189'
# LaTeX_RadioNuclide = r'$\mathrm{^{189}Re}$'


# nuclide = r'Re190'
# LaTeX_RadioNuclide = r'$\mathrm{^{190}Re}$'


# nuclide = r'Re190m'
# LaTeX_RadioNuclide = r'$\mathrm{^{190m}Re}$'


# nuclide = r'Rh94'
# LaTeX_RadioNuclide = r'$\mathrm{^{94}Rh}$'


# nuclide = r'Rh95'
# LaTeX_RadioNuclide = r'$\mathrm{^{95}Rh}$'


# nuclide = r'Rh95m'
# LaTeX_RadioNuclide = r'$\mathrm{^{95m}Rh}$'


# nuclide = r'Rh96'
# LaTeX_RadioNuclide = r'$\mathrm{^{96}Rh}$'


# nuclide = r'Rh96m'
# LaTeX_RadioNuclide = r'$\mathrm{^{96m}Rh}$'


# nuclide = r'Rh97'
# LaTeX_RadioNuclide = r'$\mathrm{^{97}Rh}$'


# nuclide = r'Rh97m'
# LaTeX_RadioNuclide = r'$\mathrm{^{97m}Rh}$'


# nuclide = r'Rh98'
# LaTeX_RadioNuclide = r'$\mathrm{^{98}Rh}$'


# nuclide = r'Rh99'
# LaTeX_RadioNuclide = r'$\mathrm{^{99}Rh}$'


# nuclide = r'Rh99m'
# LaTeX_RadioNuclide = r'$\mathrm{^{99m}Rh}$'


# nuclide = r'Rh100'
# LaTeX_RadioNuclide = r'$\mathrm{^{100}Rh}$'


# nuclide = r'Rh100m'
# LaTeX_RadioNuclide = r'$\mathrm{^{100m}Rh}$'


# nuclide = r'Rh101'
# LaTeX_RadioNuclide = r'$\mathrm{^{101}Rh}$'


# nuclide = r'Rh101m'
# LaTeX_RadioNuclide = r'$\mathrm{^{101m}Rh}$'


# nuclide = r'Rh102'
# LaTeX_RadioNuclide = r'$\mathrm{^{102}Rh}$'


# nuclide = r'Rh102m'
# LaTeX_RadioNuclide = r'$\mathrm{^{102m}Rh}$'


# nuclide = r'Rh103m'
# LaTeX_RadioNuclide = r'$\mathrm{^{103m}Rh}$'


# nuclide = r'Rh104'
# LaTeX_RadioNuclide = r'$\mathrm{^{104}Rh}$'


# nuclide = r'Rh104m'
# LaTeX_RadioNuclide = r'$\mathrm{^{104m}Rh}$'


# nuclide = r'Rh105'
# LaTeX_RadioNuclide = r'$\mathrm{^{105}Rh}$'


# nuclide = r'Rh106'
# LaTeX_RadioNuclide = r'$\mathrm{^{106}Rh}$'


# nuclide = r'Rh106m'
# LaTeX_RadioNuclide = r'$\mathrm{^{106m}Rh}$'


# nuclide = r'Rh107'
# LaTeX_RadioNuclide = r'$\mathrm{^{107}Rh}$'


# nuclide = r'Rh108'
# LaTeX_RadioNuclide = r'$\mathrm{^{108}Rh}$'


# nuclide = r'Rh109'
# LaTeX_RadioNuclide = r'$\mathrm{^{109}Rh}$'


# nuclide = r'Rn207'
# LaTeX_RadioNuclide = r'$\mathrm{^{207}Rn}$'


# nuclide = r'Rn209'
# LaTeX_RadioNuclide = r'$\mathrm{^{209}Rn}$'


# nuclide = r'Rn210'
# LaTeX_RadioNuclide = r'$\mathrm{^{210}Rn}$'


# nuclide = r'Rn211'
# LaTeX_RadioNuclide = r'$\mathrm{^{211}Rn}$'


# nuclide = r'Rn212'
# LaTeX_RadioNuclide = r'$\mathrm{^{212}Rn}$'


# nuclide = r'Rn215'
# LaTeX_RadioNuclide = r'$\mathrm{^{215}Rn}$'


# nuclide = r'Rn216'
# LaTeX_RadioNuclide = r'$\mathrm{^{216}Rn}$'


# nuclide = r'Rn217'
# LaTeX_RadioNuclide = r'$\mathrm{^{217}Rn}$'


# nuclide = r'Rn218'
# LaTeX_RadioNuclide = r'$\mathrm{^{218}Rn}$'


# nuclide = r'Rn219'
# LaTeX_RadioNuclide = r'$\mathrm{^{219}Rn}$'


# nuclide = r'Rn220'
# LaTeX_RadioNuclide = r'$\mathrm{^{220}Rn}$'


# nuclide = r'Rn222'
# LaTeX_RadioNuclide = r'$\mathrm{^{222}Rn}$'


# nuclide = r'Rn222_FE'
# LaTeX_RadioNuclide = r'$\mathrm{^{222}Rn[FE]}$'


# nuclide = r'Rn222_RE'
# LaTeX_RadioNuclide = r'$\mathrm{^{222}Rn_{RE}}$'


# nuclide = r'Rn223'
# LaTeX_RadioNuclide = r'$\mathrm{^{223}Rn}$'


# nuclide = r'Ru92'
# LaTeX_RadioNuclide = r'$\mathrm{^{92}Ru}$'


# nuclide = r'Ru94'
# LaTeX_RadioNuclide = r'$\mathrm{^{94}Ru}$'


# nuclide = r'Ru95'
# LaTeX_RadioNuclide = r'$\mathrm{^{95}Ru}$'


# nuclide = r'Ru97'
# LaTeX_RadioNuclide = r'$\mathrm{^{97}Ru}$'


# nuclide = r'Ru103'
# LaTeX_RadioNuclide = r'$\mathrm{^{103}Ru}$'


# nuclide = r'Ru105'
# LaTeX_RadioNuclide = r'$\mathrm{^{105}Ru}$'


# nuclide = r'Ru106'
# LaTeX_RadioNuclide = r'$\mathrm{^{106}Ru}$'


# nuclide = r'Ru107'
# LaTeX_RadioNuclide = r'$\mathrm{^{107}Ru}$'


# nuclide = r'Ru108'
# LaTeX_RadioNuclide = r'$\mathrm{^{108}Ru}$'


# nuclide = r'S35'
# LaTeX_RadioNuclide = r'$\mathrm{_{16}^{35}S}$'


# nuclide = r'S37'
# LaTeX_RadioNuclide = r'$\mathrm{^{37}S}$'


# nuclide = r'S38'
# LaTeX_RadioNuclide = r'$\mathrm{^{38}S}$'


# nuclide = r'Sb111'
# LaTeX_RadioNuclide = r'$\mathrm{^{111}Sb}$'


# nuclide = r'Sb113'
# LaTeX_RadioNuclide = r'$\mathrm{^{113}Sb}$'


# nuclide = r'Sb114'
# LaTeX_RadioNuclide = r'$\mathrm{^{114}Sb}$'


# nuclide = r'Sb115'
# LaTeX_RadioNuclide = r'$\mathrm{^{115}Sb}$'


# nuclide = r'Sb116'
# LaTeX_RadioNuclide = r'$\mathrm{^{116}Sb}$'


# nuclide = r'Sb116m'
# LaTeX_RadioNuclide = r'$\mathrm{^{116m}Sb}$'


# nuclide = r'Sb117'
# LaTeX_RadioNuclide = r'$\mathrm{^{117}Sb}$'


# nuclide = r'Sb118'
# LaTeX_RadioNuclide = r'$\mathrm{^{118}Sb}$'


# nuclide = r'Sb118m'
# LaTeX_RadioNuclide = r'$\mathrm{^{118m}Sb}$'


# nuclide = r'Sb119'
# LaTeX_RadioNuclide = r'$\mathrm{^{119}Sb}$'


# nuclide = r'Sb120'
# LaTeX_RadioNuclide = r'$\mathrm{^{120}Sb}$'


# nuclide = r'Sb120m'
# LaTeX_RadioNuclide = r'$\mathrm{^{120m}Sb}$'


# nuclide = r'Sb122'
# LaTeX_RadioNuclide = r'$\mathrm{^{122}Sb}$'


# nuclide = r'Sb122m'
# LaTeX_RadioNuclide = r'$\mathrm{^{122m}Sb}$'


# nuclide = r'Sb124'
# LaTeX_RadioNuclide = r'$\mathrm{^{124}Sb}$'


# nuclide = r'Sb124m'
# LaTeX_RadioNuclide = r'$\mathrm{^{124m}Sb}$'


# nuclide = r'Sb124m2'
# LaTeX_RadioNuclide = r'$\mathrm{^{124m2}Sb}$'


# nuclide = r'Sb125'
# LaTeX_RadioNuclide = r'$\mathrm{^{125}Sb}$'


# nuclide = r'Sb126'
# LaTeX_RadioNuclide = r'$\mathrm{^{126}Sb}$'


# nuclide = r'Sb126m'
# LaTeX_RadioNuclide = r'$\mathrm{^{126m}Sb}$'


# nuclide = r'Sb127'
# LaTeX_RadioNuclide = r'$\mathrm{^{127}Sb}$'


# nuclide = r'Sb128'
# LaTeX_RadioNuclide = r'$\mathrm{^{128}Sb}$'


# nuclide = r'Sb128m'
# LaTeX_RadioNuclide = r'$\mathrm{^{128m}Sb}$'


# nuclide = r'Sb129'
# LaTeX_RadioNuclide = r'$\mathrm{^{129}Sb}$'


# nuclide = r'Sb130'
# LaTeX_RadioNuclide = r'$\mathrm{^{130}Sb}$'


# nuclide = r'Sb130m'
# LaTeX_RadioNuclide = r'$\mathrm{^{130m}Sb}$'


# nuclide = r'Sb131'
# LaTeX_RadioNuclide = r'$\mathrm{^{131}Sb}$'


# nuclide = r'Sb133'
# LaTeX_RadioNuclide = r'$\mathrm{^{133}Sb}$'


# nuclide = r'Sc42m'
# LaTeX_RadioNuclide = r'$\mathrm{^{42m}Sc}$'


# nuclide = r'Sc43'
# LaTeX_RadioNuclide = r'$\mathrm{_{21}^{43}Sc}$'


# nuclide = r'Sc44'
# LaTeX_RadioNuclide = r'$\mathrm{_{21}^{44}Sc}$'


# nuclide = r'Sc44m'
# LaTeX_RadioNuclide = r'$\mathrm{^{44m}Sc}$'


# nuclide = r'Sc46'
# LaTeX_RadioNuclide = r'$\mathrm{_{21}^{46}Sc}$'


# nuclide = r'Sc47'
# LaTeX_RadioNuclide = r'$\mathrm{_{21}^{47}Sc}$'


# nuclide = r'Sc48'
# LaTeX_RadioNuclide = r'$\mathrm{^{48}Sc}$'


# nuclide = r'Sc49'
# LaTeX_RadioNuclide = r'$\mathrm{^{49}Sc}$'


# nuclide = r'Sc50'
# LaTeX_RadioNuclide = r'$\mathrm{^{50}Sc}$'


# nuclide = r'Se70'
# LaTeX_RadioNuclide = r'$\mathrm{^{70}Se}$'


# nuclide = r'Se71'
# LaTeX_RadioNuclide = r'$\mathrm{^{71}Se}$'


# nuclide = r'Se72'
# LaTeX_RadioNuclide = r'$\mathrm{^{72}Se}$'


# nuclide = r'Se73'
# LaTeX_RadioNuclide = r'$\mathrm{^{73}Se}$'


# nuclide = r'Se73m'
# LaTeX_RadioNuclide = r'$\mathrm{^{73m}Se}$'


# nuclide = r'Se75'
# LaTeX_RadioNuclide = r'$\mathrm{_{34}^{75}Se}$'


# nuclide = r'Se77m'
# LaTeX_RadioNuclide = r'$\mathrm{^{77m}Se}$'


# nuclide = r'Se79'
# LaTeX_RadioNuclide = r'$\mathrm{^{79}Se}$'


# nuclide = r'Se79m'
# LaTeX_RadioNuclide = r'$\mathrm{^{79m}Se}$'


# nuclide = r'Se81'
# LaTeX_RadioNuclide = r'$\mathrm{^{81}Se}$'


# nuclide = r'Se81m'
# LaTeX_RadioNuclide = r'$\mathrm{^{81m}Se}$'


# nuclide = r'Se83'
# LaTeX_RadioNuclide = r'$\mathrm{^{83}Se}$'


# nuclide = r'Se83m'
# LaTeX_RadioNuclide = r'$\mathrm{^{83m}Se}$'


# nuclide = r'Se84'
# LaTeX_RadioNuclide = r'$\mathrm{^{84}Se}$'


# nuclide = r'Si31'
# LaTeX_RadioNuclide = r'$\mathrm{^{31}Si}$'


# nuclide = r'Si32'
# LaTeX_RadioNuclide = r'$\mathrm{^{32}Si}$'


# nuclide = r'Sm139'
# LaTeX_RadioNuclide = r'$\mathrm{^{139}Sm}$'


# nuclide = r'Sm140'
# LaTeX_RadioNuclide = r'$\mathrm{^{140}Sm}$'


# nuclide = r'Sm141'
# LaTeX_RadioNuclide = r'$\mathrm{^{141}Sm}$'


# nuclide = r'Sm141m'
# LaTeX_RadioNuclide = r'$\mathrm{^{141m}Sm}$'


# nuclide = r'Sm142'
# LaTeX_RadioNuclide = r'$\mathrm{^{142}Sm}$'


# nuclide = r'Sm143'
# LaTeX_RadioNuclide = r'$\mathrm{^{143}Sm}$'


# nuclide = r'Sm143m'
# LaTeX_RadioNuclide = r'$\mathrm{^{143m}Sm}$'


# nuclide = r'Sm145'
# LaTeX_RadioNuclide = r'$\mathrm{^{145}Sm}$'


# nuclide = r'Sm146'
# LaTeX_RadioNuclide = r'$\mathrm{^{146}Sm}$'


# nuclide = r'Sm147'
# LaTeX_RadioNuclide = r'$\mathrm{^{147}Sm}$'


# nuclide = r'Sm148'
# LaTeX_RadioNuclide = r'$\mathrm{^{148}Sm}$'


# nuclide = r'Sm151'
# LaTeX_RadioNuclide = r'$\mathrm{^{151}Sm}$'


# nuclide = r'Sm153'
# LaTeX_RadioNuclide = r'$\mathrm{^{153}Sm}$'


# nuclide = r'Sm155'
# LaTeX_RadioNuclide = r'$\mathrm{^{155}Sm}$'


# nuclide = r'Sm156'
# LaTeX_RadioNuclide = r'$\mathrm{^{156}Sm}$'


# nuclide = r'Sm157'
# LaTeX_RadioNuclide = r'$\mathrm{^{157}Sm}$'


# nuclide = r'Sn106'
# LaTeX_RadioNuclide = r'$\mathrm{^{106}Sn}$'


# nuclide = r'Sn108'
# LaTeX_RadioNuclide = r'$\mathrm{^{108}Sn}$'


# nuclide = r'Sn109'
# LaTeX_RadioNuclide = r'$\mathrm{^{109}Sn}$'


# nuclide = r'Sn110'
# LaTeX_RadioNuclide = r'$\mathrm{^{110}Sn}$'


# nuclide = r'Sn111'
# LaTeX_RadioNuclide = r'$\mathrm{^{111}Sn}$'


# nuclide = r'Sn113'
# LaTeX_RadioNuclide = r'$\mathrm{^{113}Sn}$'


# nuclide = r'Sn113m'
# LaTeX_RadioNuclide = r'$\mathrm{^{113m}Sn}$'


# nuclide = r'Sn117m'
# LaTeX_RadioNuclide = r'$\mathrm{^{117m}Sn}$'


# nuclide = r'Sn119m'
# LaTeX_RadioNuclide = r'$\mathrm{^{119m}Sn}$'


# nuclide = r'Sn121'
# LaTeX_RadioNuclide = r'$\mathrm{^{121}Sn}$'


# nuclide = r'Sn121m'
# LaTeX_RadioNuclide = r'$\mathrm{^{121m}Sn}$'


# nuclide = r'Sn123'
# LaTeX_RadioNuclide = r'$\mathrm{^{123}Sn}$'


# nuclide = r'Sn123m'
# LaTeX_RadioNuclide = r'$\mathrm{^{123m}Sn}$'


# nuclide = r'Sn125'
# LaTeX_RadioNuclide = r'$\mathrm{^{125}Sn}$'


# nuclide = r'Sn125m'
# LaTeX_RadioNuclide = r'$\mathrm{^{125m}Sn}$'


# nuclide = r'Sn126'
# LaTeX_RadioNuclide = r'$\mathrm{^{126}Sn}$'


# nuclide = r'Sn127'
# LaTeX_RadioNuclide = r'$\mathrm{^{127}Sn}$'


# nuclide = r'Sn127m'
# LaTeX_RadioNuclide = r'$\mathrm{^{127m}Sn}$'


# nuclide = r'Sn128'
# LaTeX_RadioNuclide = r'$\mathrm{^{128}Sn}$'


# nuclide = r'Sn129'
# LaTeX_RadioNuclide = r'$\mathrm{^{129}Sn}$'


# nuclide = r'Sn130'
# LaTeX_RadioNuclide = r'$\mathrm{^{130}Sn}$'


# nuclide = r'Sn130m'
# LaTeX_RadioNuclide = r'$\mathrm{^{130m}Sn}$'


# nuclide = r'Sr79'
# LaTeX_RadioNuclide = r'$\mathrm{^{79}Sr}$'


# nuclide = r'Sr80'
# LaTeX_RadioNuclide = r'$\mathrm{^{80}Sr}$'


# nuclide = r'Sr81'
# LaTeX_RadioNuclide = r'$\mathrm{^{81}Sr}$'


# nuclide = r'Sr82'
# LaTeX_RadioNuclide = r'$\mathrm{_{38}^{82}Sr}$'


# nuclide = r'Sr83'
# LaTeX_RadioNuclide = r'$\mathrm{^{83}Sr}$'


# nuclide = r'Sr85'
# LaTeX_RadioNuclide = r'$\mathrm{_{38}^{85}Sr}$'


# nuclide = r'Sr85m'
# LaTeX_RadioNuclide = r'$\mathrm{^{85m}Sr}$'


# nuclide = r'Sr87m'
# LaTeX_RadioNuclide = r'$\mathrm{^{87m}Sr}$'


# nuclide = r'Sr89'
# LaTeX_RadioNuclide = r'$\mathrm{_{38}^{89}Sr}$'


# nuclide = r'Sr90'
# LaTeX_RadioNuclide = r'$\mathrm{^{90}Sr}$'


# nuclide = r'Sr90_Y90'
# LaTeX_RadioNuclide = r'$\mathrm{_{38}^{90}Sr} + \mathrm{_{39}^{90}Y}$'


# nuclide = r'Sr91'
# LaTeX_RadioNuclide = r'$\mathrm{^{91}Sr}$'


# nuclide = r'Sr92'
# LaTeX_RadioNuclide = r'$\mathrm{^{92}Sr}$'


# nuclide = r'Sr93'
# LaTeX_RadioNuclide = r'$\mathrm{^{93}Sr}$'


# nuclide = r'Sr94'
# LaTeX_RadioNuclide = r'$\mathrm{^{94}Sr}$'


# nuclide = r'Ta170'
# LaTeX_RadioNuclide = r'$\mathrm{^{170}Ta}$'


# nuclide = r'Ta172'
# LaTeX_RadioNuclide = r'$\mathrm{^{172}Ta}$'


# nuclide = r'Ta173'
# LaTeX_RadioNuclide = r'$\mathrm{^{173}Ta}$'


# nuclide = r'Ta174'
# LaTeX_RadioNuclide = r'$\mathrm{^{174}Ta}$'


# nuclide = r'Ta175'
# LaTeX_RadioNuclide = r'$\mathrm{^{175}Ta}$'


# nuclide = r'Ta176'
# LaTeX_RadioNuclide = r'$\mathrm{^{176}Ta}$'


# nuclide = r'Ta177'
# LaTeX_RadioNuclide = r'$\mathrm{^{177}Ta}$'


# nuclide = r'Ta178'
# LaTeX_RadioNuclide = r'$\mathrm{^{178}Ta}$'


# nuclide = r'Ta178m'
# LaTeX_RadioNuclide = r'$\mathrm{^{178m}Ta}$'


# nuclide = r'Ta179'
# LaTeX_RadioNuclide = r'$\mathrm{^{179}Ta}$'


# nuclide = r'Ta180'
# LaTeX_RadioNuclide = r'$\mathrm{^{180}Ta}$'


# nuclide = r'Ta182'
# LaTeX_RadioNuclide = r'$\mathrm{^{182}Ta}$'


# nuclide = r'Ta182m'
# LaTeX_RadioNuclide = r'$\mathrm{^{182m}Ta}$'


# nuclide = r'Ta183'
# LaTeX_RadioNuclide = r'$\mathrm{^{183}Ta}$'


# nuclide = r'Ta184'
# LaTeX_RadioNuclide = r'$\mathrm{^{184}Ta}$'


# nuclide = r'Ta185'
# LaTeX_RadioNuclide = r'$\mathrm{^{185}Ta}$'


# nuclide = r'Ta186'
# LaTeX_RadioNuclide = r'$\mathrm{^{186}Ta}$'


# nuclide = r'Tb146'
# LaTeX_RadioNuclide = r'$\mathrm{^{146}Tb}$'


# nuclide = r'Tb147'
# LaTeX_RadioNuclide = r'$\mathrm{^{147}Tb}$'


# nuclide = r'Tb147m'
# LaTeX_RadioNuclide = r'$\mathrm{^{147m}Tb}$'


# nuclide = r'Tb148'
# LaTeX_RadioNuclide = r'$\mathrm{^{148}Tb}$'


# nuclide = r'Tb148m'
# LaTeX_RadioNuclide = r'$\mathrm{^{148m}Tb}$'


# nuclide = r'Tb149'
# LaTeX_RadioNuclide = r'$\mathrm{^{149}Tb}$'


# nuclide = r'Tb149m'
# LaTeX_RadioNuclide = r'$\mathrm{^{149m}Tb}$'


# nuclide = r'Tb150'
# LaTeX_RadioNuclide = r'$\mathrm{^{150}Tb}$'


# nuclide = r'Tb150m'
# LaTeX_RadioNuclide = r'$\mathrm{^{150m}Tb}$'


# nuclide = r'Tb151'
# LaTeX_RadioNuclide = r'$\mathrm{^{151}Tb}$'


# nuclide = r'Tb151m'
# LaTeX_RadioNuclide = r'$\mathrm{^{151m}Tb}$'


# nuclide = r'Tb152'
# LaTeX_RadioNuclide = r'$\mathrm{^{152}Tb}$'


# nuclide = r'Tb152m'
# LaTeX_RadioNuclide = r'$\mathrm{^{152m}Tb}$'


# nuclide = r'Tb153'
# LaTeX_RadioNuclide = r'$\mathrm{^{153}Tb}$'


# nuclide = r'Tb154'
# LaTeX_RadioNuclide = r'$\mathrm{^{154}Tb}$'


# nuclide = r'Tb155'
# LaTeX_RadioNuclide = r'$\mathrm{^{155}Tb}$'


# nuclide = r'Tb156'
# LaTeX_RadioNuclide = r'$\mathrm{^{156}Tb}$'


# nuclide = r'Tb156m'
# LaTeX_RadioNuclide = r'$\mathrm{^{156m}Tb}$'


# nuclide = r'Tb156m2'
# LaTeX_RadioNuclide = r'$\mathrm{^{156m2}Tb}$'


# nuclide = r'Tb157'
# LaTeX_RadioNuclide = r'$\mathrm{^{157}Tb}$'


# nuclide = r'Tb158'
# LaTeX_RadioNuclide = r'$\mathrm{^{158}Tb}$'


# nuclide = r'Tb160'
# LaTeX_RadioNuclide = r'$\mathrm{^{160}Tb}$'


# nuclide = r'Tb161'
# LaTeX_RadioNuclide = r'$\mathrm{^{161}Tb}$'


# nuclide = r'Tb162'
# LaTeX_RadioNuclide = r'$\mathrm{^{162}Tb}$'


# nuclide = r'Tb163'
# LaTeX_RadioNuclide = r'$\mathrm{^{163}Tb}$'


# nuclide = r'Tb164'
# LaTeX_RadioNuclide = r'$\mathrm{^{164}Tb}$'


# nuclide = r'Tb165'
# LaTeX_RadioNuclide = r'$\mathrm{^{165}Tb}$'


# nuclide = r'Tc91'
# LaTeX_RadioNuclide = r'$\mathrm{^{91}Tc}$'


# nuclide = r'Tc91m'
# LaTeX_RadioNuclide = r'$\mathrm{^{91m}Tc}$'


# nuclide = r'Tc92'
# LaTeX_RadioNuclide = r'$\mathrm{^{92}Tc}$'


# nuclide = r'Tc93'
# LaTeX_RadioNuclide = r'$\mathrm{^{93}Tc}$'


# nuclide = r'Tc93m'
# LaTeX_RadioNuclide = r'$\mathrm{^{93m}Tc}$'


# nuclide = r'Tc94'
# LaTeX_RadioNuclide = r'$\mathrm{^{94}Tc}$'


# nuclide = r'Tc94m'
# LaTeX_RadioNuclide = r'$\mathrm{^{94m}Tc}$'


# nuclide = r'Tc95'
# LaTeX_RadioNuclide = r'$\mathrm{^{95}Tc}$'


# nuclide = r'Tc95m'
# LaTeX_RadioNuclide = r'$\mathrm{^{95m}Tc}$'


# nuclide = r'Tc96'
# LaTeX_RadioNuclide = r'$\mathrm{^{96}Tc}$'


# nuclide = r'Tc96m'
# LaTeX_RadioNuclide = r'$\mathrm{^{96m}Tc}$'


# nuclide = r'Tc97'
# LaTeX_RadioNuclide = r'$\mathrm{^{97}Tc}$'


# nuclide = r'Tc97m'
# LaTeX_RadioNuclide = r'$\mathrm{^{97m}Tc}$'


# nuclide = r'Tc98'
# LaTeX_RadioNuclide = r'$\mathrm{^{98}Tc}$'


# nuclide = r'Tc99'
# LaTeX_RadioNuclide = r'$\mathrm{_{43}^{99}Tc}$'


# nuclide = r'Tc99m'
# LaTeX_RadioNuclide = r'$\mathrm{^{99m}Tc}$'


# nuclide = r'Tc101'
# LaTeX_RadioNuclide = r'$\mathrm{^{101}Tc}$'


# nuclide = r'Tc102'
# LaTeX_RadioNuclide = r'$\mathrm{^{102}Tc}$'


# nuclide = r'Tc102m'
# LaTeX_RadioNuclide = r'$\mathrm{^{102m}Tc}$'


# nuclide = r'Tc104'
# LaTeX_RadioNuclide = r'$\mathrm{^{104}Tc}$'


# nuclide = r'Tc105'
# LaTeX_RadioNuclide = r'$\mathrm{^{105}Tc}$'


# nuclide = r'Te113'
# LaTeX_RadioNuclide = r'$\mathrm{^{113}Te}$'


# nuclide = r'Te114'
# LaTeX_RadioNuclide = r'$\mathrm{^{114}Te}$'


# nuclide = r'Te115'
# LaTeX_RadioNuclide = r'$\mathrm{^{115}Te}$'


# nuclide = r'Te115m'
# LaTeX_RadioNuclide = r'$\mathrm{^{115m}Te}$'


# nuclide = r'Te116'
# LaTeX_RadioNuclide = r'$\mathrm{^{116}Te}$'


# nuclide = r'Te117'
# LaTeX_RadioNuclide = r'$\mathrm{^{117}Te}$'


# nuclide = r'Te118'
# LaTeX_RadioNuclide = r'$\mathrm{^{118}Te}$'


# nuclide = r'Te119'
# LaTeX_RadioNuclide = r'$\mathrm{^{119}Te}$'


# nuclide = r'Te119m'
# LaTeX_RadioNuclide = r'$\mathrm{^{119m}Te}$'


# nuclide = r'Te121'
# LaTeX_RadioNuclide = r'$\mathrm{^{121}Te}$'


# nuclide = r'Te121m'
# LaTeX_RadioNuclide = r'$\mathrm{^{121m}Te}$'


# nuclide = r'Te123'
# LaTeX_RadioNuclide = r'$\mathrm{^{123}Te}$'


# nuclide = r'Te123m'
# LaTeX_RadioNuclide = r'$\mathrm{^{123m}Te}$'


# nuclide = r'Te125m'
# LaTeX_RadioNuclide = r'$\mathrm{^{125m}Te}$'


# nuclide = r'Te127'
# LaTeX_RadioNuclide = r'$\mathrm{^{127}Te}$'


# nuclide = r'Te127m'
# LaTeX_RadioNuclide = r'$\mathrm{^{127m}Te}$'


# nuclide = r'Te129'
# LaTeX_RadioNuclide = r'$\mathrm{^{129}Te}$'


# nuclide = r'Te129m'
# LaTeX_RadioNuclide = r'$\mathrm{^{129m}Te}$'


# nuclide = r'Te131'
# LaTeX_RadioNuclide = r'$\mathrm{^{131}Te}$'


# nuclide = r'Te131m'
# LaTeX_RadioNuclide = r'$\mathrm{^{131m}Te}$'


# nuclide = r'Te132'
# LaTeX_RadioNuclide = r'$\mathrm{^{132}Te}$'


# nuclide = r'Te133'
# LaTeX_RadioNuclide = r'$\mathrm{^{133}Te}$'


# nuclide = r'Te133m'
# LaTeX_RadioNuclide = r'$\mathrm{^{133m}Te}$'


# nuclide = r'Te134'
# LaTeX_RadioNuclide = r'$\mathrm{^{134}Te}$'


# nuclide = r'Th223'
# LaTeX_RadioNuclide = r'$\mathrm{^{223}Th}$'


# nuclide = r'Th224'
# LaTeX_RadioNuclide = r'$\mathrm{^{224}Th}$'


# nuclide = r'Th226'
# LaTeX_RadioNuclide = r'$\mathrm{^{226}Th}$'


# nuclide = r'Th227'
# LaTeX_RadioNuclide = r'$\mathrm{^{227}Th}$'


# nuclide = r'Th228'
# LaTeX_RadioNuclide = r'$\mathrm{^{228}Th}$'


# nuclide = r'Th229'
# LaTeX_RadioNuclide = r'$\mathrm{^{229}Th}$'


# nuclide = r'Th230'
# LaTeX_RadioNuclide = r'$\mathrm{^{230}Th}$'


# nuclide = r'Th231'
# LaTeX_RadioNuclide = r'$\mathrm{^{231}Th}$'


# nuclide = r'Th232'
# LaTeX_RadioNuclide = r'$\mathrm{^{232}Th}$'


# nuclide = r'Th232_FE'
# LaTeX_RadioNuclide = r'$\mathrm{^{232}Th[FE]}$'


# nuclide = r'Th233'
# LaTeX_RadioNuclide = r'$\mathrm{^{233}Th}$'


# nuclide = r'Th234'
# LaTeX_RadioNuclide = r'$\mathrm{^{234}Th}$'


# nuclide = r'Th235'
# LaTeX_RadioNuclide = r'$\mathrm{^{235}Th}$'


# nuclide = r'Th236'
# LaTeX_RadioNuclide = r'$\mathrm{^{236}Th}$'


# nuclide = r'Ti44'
# LaTeX_RadioNuclide = r'$\mathrm{^{44}Ti}$'


# nuclide = r'Ti45'
# LaTeX_RadioNuclide = r'$\mathrm{_{22}^{45}Sc}$'


# nuclide = r'Ti51'
# LaTeX_RadioNuclide = r'$\mathrm{^{51}Ti}$'


# nuclide = r'Ti52'
# LaTeX_RadioNuclide = r'$\mathrm{^{52}Ti}$'


# nuclide = r'Tl190'
# LaTeX_RadioNuclide = r'$\mathrm{^{190}Tl}$'


# nuclide = r'Tl190m'
# LaTeX_RadioNuclide = r'$\mathrm{^{190m}Tl}$'


# nuclide = r'Tl194'
# LaTeX_RadioNuclide = r'$\mathrm{^{194}Tl}$'


# nuclide = r'Tl194m'
# LaTeX_RadioNuclide = r'$\mathrm{^{194m}Tl}$'


# nuclide = r'Tl195'
# LaTeX_RadioNuclide = r'$\mathrm{^{195}Tl}$'


# nuclide = r'Tl196'
# LaTeX_RadioNuclide = r'$\mathrm{^{196}Tl}$'


# nuclide = r'Tl197'
# LaTeX_RadioNuclide = r'$\mathrm{^{197}Tl}$'


# nuclide = r'Tl198'
# LaTeX_RadioNuclide = r'$\mathrm{^{198}Tl}$'


# nuclide = r'Tl198m'
# LaTeX_RadioNuclide = r'$\mathrm{^{198m}Tl}$'


# nuclide = r'Tl199'
# LaTeX_RadioNuclide = r'$\mathrm{^{199}Tl}$'


# nuclide = r'Tl200'
# LaTeX_RadioNuclide = r'$\mathrm{^{200}Tl}$'


# nuclide = r'Tl201'
# LaTeX_RadioNuclide = r'$\mathrm{^{201}Tl}$'


# nuclide = r'Tl202'
# LaTeX_RadioNuclide = r'$\mathrm{^{202}Tl}$'


# nuclide = r'Tl204'
# LaTeX_RadioNuclide = r'$\mathrm{^{204}Tl}$'


# nuclide = r'Tl206'
# LaTeX_RadioNuclide = r'$\mathrm{^{206}Tl}$'


# nuclide = r'Tl206m'
# LaTeX_RadioNuclide = r'$\mathrm{^{206m}Tl}$'


# nuclide = r'Tl207'
# LaTeX_RadioNuclide = r'$\mathrm{^{207}Tl}$'


# nuclide = r'Tl208'
# LaTeX_RadioNuclide = r'$\mathrm{^{208}Tl}$'


# nuclide = r'Tl209'
# LaTeX_RadioNuclide = r'$\mathrm{^{209}Tl}$'


# nuclide = r'Tl210'
# LaTeX_RadioNuclide = r'$\mathrm{^{210}Tl}$'


# nuclide = r'Tm161'
# LaTeX_RadioNuclide = r'$\mathrm{^{161}Tm}$'


# nuclide = r'Tm162'
# LaTeX_RadioNuclide = r'$\mathrm{^{162}Tm}$'


# nuclide = r'Tm163'
# LaTeX_RadioNuclide = r'$\mathrm{^{163}Tm}$'


# nuclide = r'Tm164'
# LaTeX_RadioNuclide = r'$\mathrm{^{164}Tm}$'


# nuclide = r'Tm165'
# LaTeX_RadioNuclide = r'$\mathrm{^{165}Tm}$'


# nuclide = r'Tm166'
# LaTeX_RadioNuclide = r'$\mathrm{^{166}Tm}$'


# nuclide = r'Tm167'
# LaTeX_RadioNuclide = r'$\mathrm{^{167}Tm}$'


# nuclide = r'Tm168'
# LaTeX_RadioNuclide = r'$\mathrm{^{168}Tm}$'


# nuclide = r'Tm170'
# LaTeX_RadioNuclide = r'$\mathrm{^{170}Tm}$'


# nuclide = r'Tm171'
# LaTeX_RadioNuclide = r'$\mathrm{^{171}Tm}$'


# nuclide = r'Tm172'
# LaTeX_RadioNuclide = r'$\mathrm{^{172}Tm}$'


# nuclide = r'Tm173'
# LaTeX_RadioNuclide = r'$\mathrm{^{173}Tm}$'


# nuclide = r'Tm174'
# LaTeX_RadioNuclide = r'$\mathrm{^{174}Tm}$'


# nuclide = r'Tm175'
# LaTeX_RadioNuclide = r'$\mathrm{^{175}Tm}$'


# nuclide = r'Tm176'
# LaTeX_RadioNuclide = r'$\mathrm{^{176}Tm}$'


# nuclide = r'TuccAv_KThU_FE'
# LaTeX_RadioNuclide = r'$\text{K-Th-Unat[FE]}$'


# nuclide = r'U227'
# LaTeX_RadioNuclide = r'$\mathrm{^{227}U}$'


# nuclide = r'U228'
# LaTeX_RadioNuclide = r'$\mathrm{^{228}U}$'


# nuclide = r'U230'
# LaTeX_RadioNuclide = r'$\mathrm{^{230}U}$'


# nuclide = r'U231'
# LaTeX_RadioNuclide = r'$\mathrm{^{231}U}$'


# nuclide = r'U232'
# LaTeX_RadioNuclide = r'$\mathrm{^{232}U}$'


# nuclide = r'U233'
# LaTeX_RadioNuclide = r'$\mathrm{^{233}U}$'


# nuclide = r'U234'
# LaTeX_RadioNuclide = r'$\mathrm{^{234}U}$'


# nuclide = r'U235'
# LaTeX_RadioNuclide = r'$\mathrm{^{235}U}$'


# nuclide = r'U235m'
# LaTeX_RadioNuclide = r'$\mathrm{^{235m}U}$'


# nuclide = r'U236'
# LaTeX_RadioNuclide = r'$\mathrm{^{236}U}$'


# nuclide = r'U237'
# LaTeX_RadioNuclide = r'$\mathrm{^{237}U}$'


# nuclide = r'U238'
# LaTeX_RadioNuclide = r'$\mathrm{^{238}U}$'


# nuclide = r'U239'
# LaTeX_RadioNuclide = r'$\mathrm{^{239}U}$'


# nuclide = r'U240'
# LaTeX_RadioNuclide = r'$\mathrm{^{240}U}$'


# nuclide = r'U232_FE'
# LaTeX_RadioNuclide = r'$\mathrm{^{232}U[FE]}$'


# nuclide = r'U238_FE'
# LaTeX_RadioNuclide = r'$\mathrm{^{238}U[FE]}$'


# nuclide = r'Unat'
# LaTeX_RadioNuclide = r'$\mathrm{^{nat}U}$'


# nuclide = r'Unat_FE'
# LaTeX_RadioNuclide = r'$\mathrm{^{nat}U[FE]}$'


# nuclide = r'Urec_LEU'
# LaTeX_RadioNuclide = 'recU from LEU'


# nuclide = r'Urec_Th232_U233_FC'
# LaTeX_RadioNuclide = r'Urec from $\mathrm{^{232}Th} \rightarrow \mathrm{^{233}U}$'


# nuclide = r'U242'
# LaTeX_RadioNuclide = r'$\mathrm{^{242}U}$'


# nuclide = r'V47'
# LaTeX_RadioNuclide = r'$\mathrm{^{47}V}$'


# nuclide = r'V48'
# LaTeX_RadioNuclide = r'$\mathrm{^{48}V}$'


# nuclide = r'V49'
# LaTeX_RadioNuclide = r'$\mathrm{_23^{49}V}$'


# nuclide = r'V50'
# LaTeX_RadioNuclide = r'$\mathrm{^{50}V}$'


# nuclide = r'V52'
# LaTeX_RadioNuclide = r'$\mathrm{_{23}^{52}V}$'


# nuclide = r'V53'
# LaTeX_RadioNuclide = r'$\mathrm{^{53}V}$'


# nuclide = r'W177'
# LaTeX_RadioNuclide = r'$\mathrm{^{177}W}$'


# nuclide = r'W178'
# LaTeX_RadioNuclide = r'$\mathrm{^{178}W}$'


# nuclide = r'W179'
# LaTeX_RadioNuclide = r'$\mathrm{^{179}W}$'


# nuclide = r'W179m'
# LaTeX_RadioNuclide = r'$\mathrm{^{179m}W}$'


# nuclide = r'W181'
# LaTeX_RadioNuclide = r'$\mathrm{^{181}W}$'


# nuclide = r'W185'
# LaTeX_RadioNuclide = r'$\mathrm{^{185}W}$'


# nuclide = r'W185m'
# LaTeX_RadioNuclide = r'$\mathrm{^{185m}W}$'


# nuclide = r'W187'
# LaTeX_RadioNuclide = r'$\mathrm{^{187}W}$'


# nuclide = r'W188'
# LaTeX_RadioNuclide = r'$\mathrm{^{188}W}$'


# nuclide = r'W190'
# LaTeX_RadioNuclide = r'$\mathrm{^{190}W}$'


# nuclide = r'Xe120'
# LaTeX_RadioNuclide = r'$\mathrm{^{120}Xe}$'


# nuclide = r'Xe121'
# LaTeX_RadioNuclide = r'$\mathrm{^{121}Xe}$'


# nuclide = r'Xe122'
# LaTeX_RadioNuclide = r'$\mathrm{^{122}Xe}$'


# nuclide = r'Xe123'
# LaTeX_RadioNuclide = r'$\mathrm{^{123}Xe}$'


# nuclide = r'Xe125'
# LaTeX_RadioNuclide = r'$\mathrm{^{125}Xe}$'


# nuclide = r'Xe127'
# LaTeX_RadioNuclide = r'$\mathrm{^{127}Xe}$'


# nuclide = r'Xe127m'
# LaTeX_RadioNuclide = r'$\mathrm{^{127m}Xe}$'


# nuclide = r'Xe129m'
# LaTeX_RadioNuclide = r'$\mathrm{^{129m}Xe}$'


# nuclide = r'Xe131m'
# LaTeX_RadioNuclide = r'$\mathrm{^{131m}Xe}$'


# nuclide = r'Xe133'
# LaTeX_RadioNuclide = r'$\mathrm{^{133}Xe}$'


# nuclide = r'Xe133m'
# LaTeX_RadioNuclide = r'$\mathrm{^{133m}Xe}$'


# nuclide = r'Xe135'
# LaTeX_RadioNuclide = r'$\mathrm{^{135}Xe}$'


# nuclide = r'Xe135m'
# LaTeX_RadioNuclide = r'$\mathrm{^{135m}Xe}$'


# nuclide = r'Xe137'
# LaTeX_RadioNuclide = r'$\mathrm{^{137}Xe}$'


# nuclide = r'Xe138'
# LaTeX_RadioNuclide = r'$\mathrm{^{138}Xe}$'


# nuclide = r'Y81'
# LaTeX_RadioNuclide = r'$\mathrm{^{81}Y}$'


# nuclide = r'Y83'
# LaTeX_RadioNuclide = r'$\mathrm{^{83}Y}$'


# nuclide = r'Y83m'
# LaTeX_RadioNuclide = r'$\mathrm{^{83m}Y}$'


# nuclide = r'Y84m'
# LaTeX_RadioNuclide = r'$\mathrm{^{84m}Y}$'


# nuclide = r'Y85'
# LaTeX_RadioNuclide = r'$\mathrm{^{85}Y}$'


# nuclide = r'Y85m'
# LaTeX_RadioNuclide = r'$\mathrm{^{85m}Y}$'


# nuclide = r'Y86'
# LaTeX_RadioNuclide = r'$\mathrm{^{86}Y}$'


# nuclide = r'Y86m'
# LaTeX_RadioNuclide = r'$\mathrm{^{86m}Y}$'


# nuclide = r'Y87'
# LaTeX_RadioNuclide = r'$\mathrm{^{87}Y}$'


# nuclide = r'Y87m'
# LaTeX_RadioNuclide = r'$\mathrm{^{87m}Y}$'


# nuclide = r'Y88'
# LaTeX_RadioNuclide = r'$\mathrm{^{88}Y}$'


# nuclide = r'Y89m'
# LaTeX_RadioNuclide = r'$\mathrm{^{89m}Y}$'


# nuclide = r'Y90'
# LaTeX_RadioNuclide = r'$\mathrm{_{39}^{90}Y}$'


# nuclide = r'Y90m'
# LaTeX_RadioNuclide = r'$\mathrm{^{90m}Y}$'


# nuclide = r'Y91'
# LaTeX_RadioNuclide = r'$\mathrm{^{91}Y}$'


# nuclide = r'Y91m'
# LaTeX_RadioNuclide = r'$\mathrm{^{91m}Y}$'


# nuclide = r'Y92'
# LaTeX_RadioNuclide = r'$\mathrm{^{92}Y}$'


# nuclide = r'Y93'
# LaTeX_RadioNuclide = r'$\mathrm{^{93}Y}$'


# nuclide = r'Y94'
# LaTeX_RadioNuclide = r'$\mathrm{^{94}Y}$'


# nuclide = r'Y95'
# LaTeX_RadioNuclide = r'$\mathrm{^{95}Y}$'


# nuclide = r'Yb162'
# LaTeX_RadioNuclide = r'$\mathrm{^{162}Yb}$'


# nuclide = r'Yb163'
# LaTeX_RadioNuclide = r'$\mathrm{^{163}Yb}$'


# nuclide = r'Yb164'
# LaTeX_RadioNuclide = r'$\mathrm{^{164}Yb}$'


# nuclide = r'Yb165'
# LaTeX_RadioNuclide = r'$\mathrm{^{165}Yb}$'


# nuclide = r'Yb166'
# LaTeX_RadioNuclide = r'$\mathrm{^{166}Yb}$'


# nuclide = r'Yb167'
# LaTeX_RadioNuclide = r'$\mathrm{^{167}Yb}$'


# nuclide = r'Yb169'
# LaTeX_RadioNuclide = r'$\mathrm{^{169}Yb}$'


# nuclide = r'Yb175'
# LaTeX_RadioNuclide = r'$\mathrm{^{175}Yb}$'


# nuclide = r'Yb177'
# LaTeX_RadioNuclide = r'$\mathrm{^{177}Yb}$'


# nuclide = r'Yb178'
# LaTeX_RadioNuclide = r'$\mathrm{^{178}Yb}$'


# nuclide = r'Yb179'
# LaTeX_RadioNuclide = r'$\mathrm{^{179}Yb}$'


# nuclide = r'Zn60'
# LaTeX_RadioNuclide = r'$\mathrm{^{60}Zn}$'


# nuclide = r'Zn61'
# LaTeX_RadioNuclide = r'$\mathrm{^{61}Zn}$'


# nuclide = r'Zn62'
# LaTeX_RadioNuclide = r'$\mathrm{^{62}Zn}$'


# nuclide = r'Zn63'
# LaTeX_RadioNuclide = r'$\mathrm{^{63}Zn}$'


# nuclide = r'Zn65'
# LaTeX_RadioNuclide = r'$\mathrm{_{30}^{65}Zn}$'


# nuclide = r'Zn69'
# LaTeX_RadioNuclide = r'$\mathrm{^{69}Zn}$'


# nuclide = r'Zn69m'
# LaTeX_RadioNuclide = r'$\mathrm{^{69m}Zn}$'


# nuclide = r'Zn71'
# LaTeX_RadioNuclide = r'$\mathrm{^{71}Zn}$'


# nuclide = r'Zn71m'
# LaTeX_RadioNuclide = r'$\mathrm{^{71m}Zn}$'


# nuclide = r'Zn72'
# LaTeX_RadioNuclide = r'$\mathrm{^{72}Zn}$'


# nuclide = r'Zr85'
# LaTeX_RadioNuclide = r'$\mathrm{^{85}Zr}$'


# nuclide = r'Zr86'
# LaTeX_RadioNuclide = r'$\mathrm{^{86}Zr}$'


# nuclide = r'Zr87'
# LaTeX_RadioNuclide = r'$\mathrm{^{87}Zr}$'


# nuclide = r'Zr88'
# LaTeX_RadioNuclide = r'$\mathrm{^{88}Zr}$'


# nuclide = r'Zr89'
# LaTeX_RadioNuclide = r'$\mathrm{^{89}Zr}$'


# nuclide = r'Zr89m'
# LaTeX_RadioNuclide = r'$\mathrm{^{89m}Zr}$'


# nuclide = r'Zr93'
# LaTeX_RadioNuclide = r'$\mathrm{_{40}^{93}Zr}$'


# nuclide = r'Zr95'
# LaTeX_RadioNuclide = r'$\mathrm{^{95}Zr}$'


# nuclide = r'Zr95_Nb95'
# LaTeX_RadioNuclide = r'$\mathrm{_{40}^{95}Zr}$ \text{and} \mathrm{_{41}^{95}Nb}$'


# nuclide = r'Zr97'
# LaTeX_RadioNuclide = r'$\mathrm{^{97}Zr}$'

#####

SpreadSheet = path_inp_abs + nuclide + r'.xlsx'

#####

# Extract relevant parameters for the specified radionuclide:
cols_req = [0,1,2,3,4,5,6,7,8,9,10]
df_Parameters = pd.read_excel(SpreadSheet,
                           sheet_name='Parameters',
                           usecols=cols_req,
                           dtype=None,
                           header=[1],
                           index_col=None,
                           nrows=None)
Name_Nuclide = df_Parameters['Nuclide'][0]
Nuclide_LaTeX = df_Parameters['Name_Nuclide_LaTeX'][0]
Activity_total = df_Parameters['Activity_tot'].to_numpy()
Energy_max = df_Parameters['Energy_max'].to_numpy()
Yield_min = df_Parameters['Yield_min'].to_numpy()
Yield_max = df_Parameters['Yield_max'].to_numpy()
Z = df_Parameters['Z'].to_numpy()
A = df_Parameters['A'].to_numpy()
ISO = df_Parameters['ISO'].to_numpy()
T_half = df_Parameters['T_half'].to_numpy()
Time_Units = df_Parameters['Time_Units'][0]
print('Activity_total  = ', Activity_total)
print('Energy_max   = ', Energy_max)
print('Yield_min   = ', Yield_min)
print('Yield_max   = ', Yield_max)


# Yields of different categories of ionising radiation emissions: SFN, DEP, DEE, BEN, BEP, DEA, DER, DEFFL, DEFFH
cols_req = [0,1,2,3,4,5,6,7,8]
df_Emission_Yields = pd.read_excel(SpreadSheet,
                                sheet_name='Yields',
                                dtype=np.float64,
                                usecols=cols_req,
                                header=[1],
                                index_col=None,
                                nrows=None)
Yield_SFN    =  df_Emission_Yields['Y_SFN'].to_numpy()
Yield_DEP    =  df_Emission_Yields['Y_DEP'].to_numpy()
Yield_DEE    =  df_Emission_Yields['Y_DEE'].to_numpy()
Yield_BEN    =  df_Emission_Yields['Y_BEN'].to_numpy()
Yield_BEP    =  df_Emission_Yields['Y_BEP'].to_numpy()
Yield_DEA    =  df_Emission_Yields['Y_DEA'].to_numpy()
Yield_DER    =  df_Emission_Yields['Y_DER'].to_numpy()
Yield_DEFFL  =  df_Emission_Yields['Y_DEFFL'].to_numpy()
Yield_DEFFH  =  df_Emission_Yields['Y_DEFFH'].to_numpy()


print('Yield_SFN    = ', Yield_SFN)
print('Yield_DEP    = ', Yield_DEP)
print('Yield_DEE    = ', Yield_DEE)
print('Yield_BEN    = ', Yield_BEN)
print('Yield_BEP    = ', Yield_BEP)
print('Yield_DEA    = ', Yield_DEA)
print('Yield_DER    = ', Yield_DER)
print('Yield_DEFFL  = ', Yield_DEFFL)
print('Yield_DEFFH  = ', Yield_DEFFH)


# Specify appropriate major tick steps:
def x_major_tick_step(Energy_max):
    if Energy_max >= 0.01:
        x_major_ts = 0.001
    if Energy_max >= 0.02:
        x_major_ts = 0.002
    if Energy_max >= 0.03:
        x_major_ts = 0.03
    if Energy_max >= 0.05:
        x_major_ts = 0.01
    if Energy_max >= 0.2:
        x_major_ts = 0.02
    if Energy_max >= 0.4:
        x_major_ts = 0.05
    if Energy_max >= 0.5:
        x_major_ts = 0.10
    if Energy_max >= 1.0:
        x_major_ts = 0.20
    if Energy_max >= 2.0:
        x_major_ts = 0.25
    if Energy_max >= 4.0:
        x_major_ts = 0.50
    if Energy_max >= 6.0:
        x_major_ts = 1.00
    if Energy_max >= 12.0:
        x_major_ts = 2.00
    return x_major_ts


print('x_major_tick_step is:', x_major_tick_step(Energy_max))


if Yield_SFN != 0.0:
    cols_req = [0, 1]
    source_spectrum_SFN = pd.read_excel(SpreadSheet,
                                        sheet_name='SFN',
                                        dtype=np.float64,
                                        usecols=cols_req,
                                        header=[1],
                                        index_col=None,
                                        nrows=None)

if Yield_DEP != 0.0:
    cols_req = [0, 1]
    source_spectrum_DEP = pd.read_excel(SpreadSheet,
                                        sheet_name='DEP',
                                        dtype=np.float64,
                                        usecols=cols_req,
                                        header=[1],
                                        index_col=None,
                                        nrows=None)

if Yield_DEE != 0.0:
    cols_req = [0, 1]
    source_spectrum_DEE = pd.read_excel(SpreadSheet,
                                        sheet_name='DEE',
                                        dtype=np.float64,
                                        usecols=cols_req,
                                        header=[1],
                                        index_col=None,
                                        nrows=None)

if Yield_BEN != 0.0:
    cols_req = [0, 1]
    source_spectrum_BEN = pd.read_excel(SpreadSheet,
                                        sheet_name='BEN',
                                        dtype=np.float64,
                                        usecols=cols_req,
                                        header=[1],
                                        index_col=None,
                                        nrows=None)

if Yield_BEP != 0.0:
    cols_req = [0, 1]
    source_spectrum_BEP = pd.read_excel(SpreadSheet,
                                        sheet_name='BEP',
                                        dtype=np.float64,
                                        usecols=cols_req,
                                        header=[1],
                                        index_col=None,
                                        nrows=None)

if Yield_DEA != 0.0:
    cols_req = [0, 1]
    source_spectrum_DEA = pd.read_excel(SpreadSheet,
                                        sheet_name='DEA',
                                        dtype=np.float64,
                                        usecols=cols_req,
                                        header=[1],
                                        index_col=None,
                                        nrows=None)



if Yield_DER != 0.0:
    cols_req = [0, 1]
    source_spectrum_DER = pd.read_excel(SpreadSheet,
                                        sheet_name='DER',
                                        dtype=np.float64,
                                        usecols=cols_req,
                                        header=[1],
                                        index_col=None,
                                        nrows=None)

if Yield_DEFFL != 0.0:
    cols_req = [0, 1]
    source_spectrum_DEFFL = pd.read_excel(SpreadSheet,
                                        sheet_name='DEFFL',
                                        dtype=np.float64,
                                        usecols=cols_req,
                                        header=[1],
                                        index_col=None,
                                        nrows=None)


if Yield_DEFFH != 0.0:
    cols_req = [0, 1]
    source_spectrum_DEFFH = pd.read_excel(SpreadSheet,
                                        sheet_name='DEFFH',
                                        dtype=np.float64,
                                        usecols=cols_req,
                                        header=[1],
                                        index_col=None,
                                        nrows=None)


if Yield_SFN != 0.0:
    eng_SFN    = source_spectrum_SFN['eng_SFN'].to_numpy()
    Q_SFN_diff = source_spectrum_SFN['pdf_yield_SFN'].to_numpy()

if Yield_DEP != 0.0:
    eng_DEP    = np.flipud(source_spectrum_DEP['eng_DEP'].to_numpy())
    Q_DEP_line = np.flipud(source_spectrum_DEP['yield_DEP'].to_numpy())

if Yield_DEE != 0.0:
    eng_DEE    = np.flipud(source_spectrum_DEE['eng_DEE'].to_numpy())
    Q_DEE_line = np.flipud(source_spectrum_DEE['yield_DEE'].to_numpy())

if Yield_BEN != 0.0:
    eng_BEN    = source_spectrum_BEN['eng_BEN'].to_numpy()
    Q_BEN_diff = source_spectrum_BEN['pdf_yield_BEN'].to_numpy()

if Yield_BEP != 0.0:
    eng_BEP    = source_spectrum_BEP['eng_BEP'].to_numpy()
    Q_BEP_diff = source_spectrum_BEP['pdf_yield_BEP'].to_numpy()

if Yield_DEA != 0.0:
    eng_DEA    = source_spectrum_DEA['eng_DEA'].to_numpy()
    Q_DEA_line = source_spectrum_DEA['yield_DEA'].to_numpy()

if Yield_DER != 0.0:
    eng_DER    = source_spectrum_DER['eng_DER'].to_numpy()
    Q_DER_line = source_spectrum_DER['yield_DER'].to_numpy()

if Yield_DEFFL != 0.0:
    eng_DEFFL    = source_spectrum_DEFFL['eng_DEFFL'].to_numpy()
    Q_DEFFL_line = source_spectrum_DEFFL['yield_DEFFL'].to_numpy()

if Yield_DEFFH != 0.0:
    eng_DEFFH    = source_spectrum_DEFFH['eng_DEFFH'].to_numpy()
    Q_DEFFH_line = source_spectrum_DEFFH['yield_DEFFH'].to_numpy()


# Correction factor: DEP
if Yield_DEP != 0.0:
    cf_DEP = (Yield_DEP / np.sum(Q_DEP_line))
    print('Correction factor for DEP = ', cf_DEP)

# Correction factor: DEE
if Yield_DEE != 0.0:
    cf_DEE = (Yield_DEE / np.sum(Q_DEE_line))
    print('Correction factor for DEE = ', cf_DEE)

# Correction factor: DEA
if Yield_DEA != 0.0:
    cf_DEA = (Yield_DEA / np.sum(Q_DEA_line))
    print('Correction factor for DEA = ', cf_DEA)

# Correction factor: DER
if Yield_DER != 0.0:
    cf_DER = (Yield_DER / np.sum(Q_DER_line))
    print('Correction factor for DER = ', cf_DER)

# Correction factor: DEFFL
if Yield_DEFFL != 0.0:
    cf_DEFFL = (Yield_DEFFL / np.sum(Q_DEFFL_line))
    print('Correction factor for DEFFL = ', cf_DEFFL)

# Correction factor: DEFFH
if Yield_DEFFH != 0.0:
    cf_DEFFH = (Yield_DEFFH / np.sum(Q_DEFFH_line))
    print('Correction factor for DEFFH = ', cf_DEFFH)


# Numerical Integration of reconstructed differential energy-spectra:
if Yield_SFN != 0.0:
    x = eng_SFN
    y = Q_SFN_diff
    if len(x) >= 2:
        spec_interp = PchipInterpolator(x, y, extrapolate=None)
        x_min = np.min(x)
        x_max = np.max(x)
        int_SFN = integrate.quad(spec_interp, x_min, x_max,
                                 epsabs=1.00E-12, epsrel=1.00E-12, limit=99999)
        int_SFN = int_SFN[0]
        print('Integral_SFN = ', int_SFN)
        cf_SFN = Yield_SFN / int_SFN
    else:
        # Single-point SFN (e.g. an average-energy line): no spectrum to
        # integrate, so the stored value already is the yield (cf = 1.0).
        cf_SFN = 1.0
    print('Correction factor for SFN = ', cf_SFN)

# Numerical Integration of reconstructed differential energy-spectra:
if Yield_BEN != 0.0:
    x = eng_BEN
    y = Q_BEN_diff
    spec_interp = PchipInterpolator(x, y, extrapolate=None)
    x_min = np.min(x)
    x_max = np.max(x)
    int_BEN = integrate.quad(spec_interp, x_min, x_max,
                             epsabs=1.00E-12, epsrel=1.00E-12, limit=99999)
    int_BEN = int_BEN[0]
    print('Integral_BEN = ', int_BEN)
    cf_BEN = Yield_BEN / int_BEN
    print('Correction factor for BEN = ', cf_BEN)

# Numerical Integration of reconstructed differential energy-spectra:
if Yield_BEP != 0.0:
    x = eng_BEP
    y = Q_BEP_diff
    spec_interp = PchipInterpolator(x, y, extrapolate=None)
    x_min = np.min(x)
    x_max = np.max(x)
    int_BEP = integrate.quad(spec_interp, x_min, x_max,
                             epsabs=1.00E-12, epsrel=1.00E-12, limit=99999)
    int_BEP = int_BEP[0]
    print('Integral_BEP = ', int_BEP)
    cf_BEP = Yield_BEP / int_BEP
    print('Correction factor for BEP = ', cf_BEP)


# Multiply Emission spectra by correction factor
if Yield_SFN != 0.0:
    Q_SFN = Q_SFN_diff * cf_SFN  # * Activity_total
if Yield_DEP != 0.0:
    Q_DEP = Q_DEP_line * cf_DEP  # * Activity_total
if Yield_DEE != 0.0:
    Q_DEE = Q_DEE_line * cf_DEE  # * Activity_total
if Yield_BEN != 0.0:
    Q_BEN = Q_BEN_diff * cf_BEN  # * Activity_total
if Yield_BEP != 0.0:
    Q_BEP = Q_BEP_diff * cf_BEP  # * Activity_total
if Yield_DEA != 0.0:
    Q_DEA = Q_DEA_line * cf_DEA  # * Activity_total
if Yield_DER != 0.0:
    Q_DER = Q_DER_line * cf_DER  # * Activity_total
if Yield_DEFFL != 0.0:
    Q_DEFFL = Q_DEFFL_line * cf_DEFFL  # * Activity_total
if Yield_DEFFH != 0.0:
    Q_DEFFH = Q_DEFFH_line * cf_DEFFH  # * Activity_total


# Plot differential SFN spectrum vs energy, as line-plot:
if Yield_SFN != 0.0:
    x = eng_SFN
    y = Q_SFN
    fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)
    fig.patch.set_linewidth(3)
    fig.patch.set_edgecolor('xkcd:black')
    axes.plot(x, y, label=r'$Q_\mathrm{{SFN}}(E)$',
              color=color_SFN, linestyle='solid', marker=None)
    #
    title  = LaTeX_RadioNuclide + r'  Spontaneous fission neutron spectrum $Q_\mathrm{{SFN}}(E)$'
    xlabel = 'Energy (MeV)'
    ylabel = r'Radiation emission spectrum $\left(\mathrm{s}^{-1}\;\mathrm{MeV}^{-1}\right)$'
    #
    axes.set_title(title, y=position_title_y)
    axes.set_xlabel(xlabel, labelpad=labelpad_x)
    axes.set_ylabel(ylabel, labelpad=labelpad_y)
    axes.set_xlim([1.0E-11, Energy_max])
    axes.set_ylim([Yield_min, Yield_max])
    axes.set_xscale(x_scale)
    axes.set_yscale(y_scale)
    axes.tick_params(axis='both', which='major')
    axes.legend(loc='best',  ncol=1, framealpha=0.95)
    axes.grid(visible=True, which='major', axis='x',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='x',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.grid(visible=True, which='major', axis='y',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='y',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.minorticks_off()
    # fig.savefig(path_out_abs + nuclide + '_SFN_line.svg', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='svg', transparent=True, pad_inches=0.00)
    # fig.savefig(path_out_abs + nuclide + '_SFN_line.pdf', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='pdf', transparent=True, pad_inches=0.00)
    # plt.show()
    plt.clf()
    plt.close()


# Plot discrete energy photon (DEP) emission spectrum vs energy as stem-plot:
if Yield_DEP != 0.0:
    x = eng_DEP
    y = Q_DEP

    def LineW_stem(x):
        if len(x) >= 0:
            lw_stem = 4.0
        if len(x) >= 50:
            lw_stem = 3.5
        if len(x) >= 250:
            lw_stem = 3.0
        if len(x) >= 500:
            lw_stem = 2.5
        if len(x) >= 750:
            lw_stem = 2.0
        if len(x) >= 1500:
            lw_stem = 1.50
        if len(x) >= 4000:
            lw_stem = 1.00
        return lw_stem

    def marker_size(x):
        if len(x) >= 0:
            marker_s = 15
        if len(x) >= 10:
            marker_s = 14
        if len(x) >= 25:
            marker_s = 13
        if len(x) >= 50:
            marker_s = 12
        if len(x) >= 250:
            marker_s = 11
        if len(x) >= 500:
            marker_s = 10
        return marker_s

    fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)
    fig.patch.set_linewidth(3)
    fig.patch.set_edgecolor('xkcd:black')
    #
    (markers, stemlines, baseline) = plt.stem(
        x, y, label=r'$Q_{\gamma}(E)$', bottom=Yield_min)
    plt.setp(markers, marker='D', markersize=marker_size(x), markerfacecolor=color_markerface_DEP,
             markeredgecolor=color_markeredge_DEP, markeredgewidth=2)
    plt.setp(stemlines, linestyle='solid', color=color_stemline_DEP,
             linewidth=LineW_stem(x), alpha=1.0)
    plt.setp(baseline, linestyle='solid', color=color_baseline,
             linewidth=lw_baseline, visible=True)
    #
    title  = LaTeX_RadioNuclide + r'  Discrete energy photon (DEP) emission spectrum $Q_{\gamma}(E)$'
    xlabel = 'Energy (MeV)'
    ylabel = 'Particle emission yield per nt'
    #
    axes.set_title(title,   y=position_title_y)
    axes.set_xlabel(xlabel, labelpad=labelpad_x)
    axes.set_ylabel(ylabel, labelpad=labelpad_y)
    axes.set_xlim([1.0E-5, Energy_max])
    axes.set_ylim([Yield_min, Yield_max])
    axes.set_xscale(x_scale)
    axes.set_yscale(y_scale)
    axes.legend(loc='best', ncol=1, framealpha=0.95)
    axes.tick_params(axis='both', which='major')
    axes.grid(visible=True, which='major', axis='x',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='x',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.grid(visible=True, which='major', axis='y',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='y',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.minorticks_off()
    # fig.savefig(path_out_abs + nuclide + '_DEP.svg', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='svg', transparent=True, pad_inches=0.00)
    # fig.savefig(path_out_abs + nuclide + '_DEP.pdf', dpi=None, facecolor='white',
    #            edgecolor='black', orientation=None, format='pdf', transparent=True, pad_inches=0.00)
    # plt.show()
    plt.clf()
    plt.close()


# Plot discrete energy electron (DEE) spectrum vs energy as stem-plot:
if Yield_DEE != 0.0:
    x = eng_DEE
    y = Q_DEE

    def LineW_stem(x):
        if len(x) >= 0:
            lw_stem = 4.0
        if len(x) >= 50:
            lw_stem = 3.5
        if len(x) >= 250:
            lw_stem = 3.0
        if len(x) >= 500:
            lw_stem = 2.5
        if len(x) >= 750:
            lw_stem = 2.0
        if len(x) >= 1500:
            lw_stem = 1.50
        if len(x) >= 4000:
            lw_stem = 1.00
        return lw_stem

    def marker_size(x):
        if len(x) >= 0:
            marker_s = 15
        if len(x) >= 10:
            marker_s = 14
        if len(x) >= 25:
            marker_s = 13
        if len(x) >= 50:
            marker_s = 12
        if len(x) >= 250:
            marker_s = 11
        if len(x) >= 500:
            marker_s = 10
        return marker_s

    fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)
    fig.patch.set_linewidth(3)
    fig.patch.set_edgecolor('xkcd:black')
    #
    (markers, stemlines, baseline) = plt.stem(
        x, y, label=r'$Q_{e}(E)$', bottom=Yield_min)
    plt.setp(markers, marker='o', markersize=marker_size(x), markerfacecolor=color_markerface_DEE,
             markeredgecolor=color_markeredge_DEE, markeredgewidth=2)
    plt.setp(stemlines, linestyle='solid', color=color_stemline_DEE,
             alpha=1.0, linewidth=LineW_stem(x))
    plt.setp(baseline, linestyle='solid', color=color_baseline,
             linewidth=lw_baseline, visible=True)
    #
    title  = LaTeX_RadioNuclide + r'  Discrete energy electron (DEE) emission spectrum $Q_{e}(E)$'
    xlabel = 'Energy (MeV)'
    ylabel = 'Particle emission yield per nt'
    #
    axes.set_title(title,   y=position_title_y)
    axes.set_xlabel(xlabel, labelpad=labelpad_x)
    axes.set_ylabel(ylabel, labelpad=labelpad_y)
    axes.set_xlim([1.0E-5, Energy_max])
    axes.set_ylim([Yield_min, Yield_max])
    axes.set_xscale(x_scale)
    axes.set_yscale(y_scale)
    axes.legend(loc='best',  ncol=1, framealpha=0.95)
    axes.tick_params(axis='both', which='major')
    axes.grid(visible=True, which='major', axis='x',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='x',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.grid(visible=True, which='major', axis='y',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='y',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.minorticks_off()
    # fig.savefig(path_out_abs + nuclide + '_DEE.svg', dpi=None, facecolor='white', edgecolor='black',
    #             orientation=None, format='svg', transparent=True, pad_inches=0.00)
    # fig.savefig(path_out_abs + nuclide + '_DEE.pdf', dpi=None, facecolor='white', edgecolor='black',
    #             orientation=None, format='pdf', transparent=True, pad_inches=0.00)
    # plt.show()
    plt.clf()
    plt.close()


# Plot differential beta-electron spectrum vs energy, as line-plot:
if Yield_BEN != 0.0:
    x = eng_BEN
    y = Q_BEN
    fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)
    fig.patch.set_linewidth(3)
    fig.patch.set_edgecolor('xkcd:black')
    axes.plot(x, y, label=r'$Q_{\beta^{-}}(E)$',
              color=color_BEN, linestyle='solid', marker=None)
    #
    title  = LaTeX_RadioNuclide + r'  Beta-electron (BEN) emission spectrum $Q_{\beta^{-}}(E)$'
    xlabel = 'Energy (MeV)'
    ylabel = r'Particle emission spectrum $\left(\mathrm{s}^{-1}\;\mathrm{MeV}^{-1}\right)$'
    #
    axes.set_title(title, y=position_title_y)
    axes.set_xlabel(xlabel, labelpad=labelpad_x)
    axes.set_ylabel(ylabel, labelpad=labelpad_y)
    axes.set_xlim([1.0E-5, Energy_max])
    axes.set_ylim([Yield_min, Yield_max])
    axes.set_xscale(x_scale)
    axes.set_yscale(y_scale)
    axes.tick_params(axis='both', which='major')
    axes.legend(loc='best',  ncol=1, framealpha=0.95)
    axes.grid(visible=True, which='major', axis='x',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='x',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.grid(visible=True, which='major', axis='y',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='y',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.minorticks_off()
    # fig.savefig(path_out_abs + nuclide + '_BEN_line.svg', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='svg', transparent=True, pad_inches=0.00)
    # fig.savefig(path_out_abs + nuclide + '_BEN_line.pdf', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='pdf', transparent=True, pad_inches=0.00)
    # plt.show()
    plt.clf()
    plt.close()


# Plot differential beta-positron spectrum vs energy, as line-plot:
if Yield_BEP != 0.0:
    x = eng_BEP
    y = Q_BEP
    fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)
    fig.patch.set_linewidth(3)
    fig.patch.set_edgecolor('xkcd:black')
    axes.plot(x, y, label=r'$Q_{\beta^{+}}(E)$',
              color=color_BEP, linestyle='solid', marker=None)
    #
    title  = LaTeX_RadioNuclide + r'  Beta-positron (BEP) emission spectrum, $Q_{\beta^{+}}(E)$'
    xlabel = 'Energy (MeV)'
    ylabel = r'Particle emission spectrum $\left(\mathrm{s}^{-1}\;\mathrm{MeV}^{-1}\right)$'
    #
    axes.set_title(title, y=position_title_y)
    axes.set_xlabel(xlabel, labelpad=labelpad_x)
    axes.set_ylabel(ylabel, labelpad=labelpad_y)
    axes.set_xlim([1.0E-5, Energy_max])
    axes.set_ylim([Yield_min, Yield_max])
    axes.set_xscale(x_scale)
    axes.set_yscale(y_scale)
    axes.tick_params(axis='both', which='major')
    axes.legend(loc='best',  ncol=1, framealpha=0.95)
    axes.grid(visible=True, which='major', axis='x',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='x',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.grid(visible=True, which='major', axis='y',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='y',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.minorticks_off()
    # fig.savefig(path_out_abs + nuclide + '_BEP_line.svg', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='svg', transparent=True, pad_inches=0.00)
    # fig.savefig(path_out_abs + nuclide + '_BEP_line.pdf', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='pdf', transparent=True, pad_inches=0.00)
    # plt.show()
    plt.clf()
    plt.close()


# Plot discrete energy alpha (DEA) spectrum vs energy as stem-plot:
if Yield_DEA != 0.0:
    x = eng_DEA
    y = Q_DEA

    def LineW_stem(x):
        if len(x) >= 0:
            lw_stem = 4.0
        if len(x) >= 50:
            lw_stem = 3.5
        if len(x) >= 250:
            lw_stem = 3.0
        if len(x) >= 500:
            lw_stem = 2.5
        if len(x) >= 750:
            lw_stem = 2.0
        if len(x) >= 1500:
            lw_stem = 1.50
        if len(x) >= 4000:
            lw_stem = 1.00
        return lw_stem

    def marker_size(x):
        if len(x) >= 0:
            marker_s = 15
        if len(x) >= 10:
            marker_s = 14
        if len(x) >= 25:
            marker_s = 13
        if len(x) >= 50:
            marker_s = 12
        if len(x) >= 250:
            marker_s = 11
        if len(x) >= 500:
            marker_s = 10
        return marker_s

    fig, axes = plt.subplots(nrows=1, ncols=1, sharex='col', sharey=True)
    (markers, stemlines, baseline) = plt.stem(
        x, y, label=r'$Q_{\alpha}(E)$', bottom=Yield_min)
    fig.patch.set_linewidth(3)
    fig.patch.set_edgecolor('xkcd:black')
    plt.setp(markers, marker='o', markersize=marker_size(x), markerfacecolor=color_markerface_DEA,
             markeredgecolor=color_markeredge_DEA, markeredgewidth=2)
    plt.setp(stemlines, linestyle='solid', color=color_stemline_DEA,
             alpha=1.0, linewidth=LineW_stem(x))
    plt.setp(baseline, linestyle='solid', color=color_baseline,
             linewidth=lw_baseline, visible=True)
    #
    title  = LaTeX_RadioNuclide + r'  Discrete energy $\alpha$-particle (DEA) emission spectrum $Q_{\alpha}(E)$'
    xlabel = 'Energy (MeV)'
    ylabel = 'Particle emission yield per nt'
    #
    axes.set_title(title, y=position_title_y)
    axes.set_xlabel(xlabel, labelpad=labelpad_x)
    axes.set_ylabel(ylabel, labelpad=labelpad_y)
    axes.set_xlim([1.0E-3, Energy_max])
    axes.set_ylim([Yield_min, Yield_max])
    axes.set_xscale(x_scale)
    axes.set_yscale(y_scale)
    axes.legend(loc='best',  ncol=1, framealpha=0.95)
    axes.tick_params(axis='both', which='major')
    axes.grid(visible=True, which='major', axis='x',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='x',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.grid(visible=True, which='major', axis='y',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='y',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.minorticks_off()
    # fig.savefig(path_out_abs + nuclide + '_DEA.svg', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='svg', transparent=True, pad_inches=0.00)
    # fig.savefig(path_out_abs + nuclide + '_DEA.pdf', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='pdf', transparent=True, pad_inches=0.00)
    # plt.show()
    plt.clf()
    plt.close()


# Plot discrete energy alpha (DER) spectrum vs energy as stem-plot:
if Yield_DER != 0.0:
    x = eng_DER
    y = Q_DER

    def LineW_stem(x):
        if len(x) >= 0:
            lw_stem = 4.0
        if len(x) >= 50:
            lw_stem = 3.5
        if len(x) >= 250:
            lw_stem = 3.0
        if len(x) >= 500:
            lw_stem = 2.5
        if len(x) >= 750:
            lw_stem = 2.0
        if len(x) >= 1500:
            lw_stem = 1.50
        if len(x) >= 4000:
            lw_stem = 1.00
        return lw_stem

    def marker_size(x):
        if len(x) >= 0:
            marker_s = 15
        if len(x) >= 10:
            marker_s = 14
        if len(x) >= 25:
            marker_s = 13
        if len(x) >= 50:
            marker_s = 12
        if len(x) >= 250:
            marker_s = 11
        if len(x) >= 500:
            marker_s = 10
        return marker_s

    fig, axes = plt.subplots(nrows=1, ncols=1, sharex='col', sharey=True)
    (markers, stemlines, baseline) = plt.stem(
        x, y, label=r'$Q_{\mathsf{rec}}(E)$', bottom=Yield_min)
    fig.patch.set_linewidth(3)
    fig.patch.set_edgecolor('xkcd:black')
    plt.setp(markers, marker='h', markersize=marker_size(x), markerfacecolor=color_markerface_DER,
             markeredgecolor=color_markeredge_DER, markeredgewidth=2)
    plt.setp(stemlines, linestyle='solid', color=color_stemline_DER,
             alpha=1.0, linewidth=LineW_stem(x))
    plt.setp(baseline, linestyle='solid', color=color_baseline,
             linewidth=lw_baseline, visible=True)
    #
    title  = LaTeX_RadioNuclide + r'  Discrete energy recoil (DER) emission spectrum $Q_{\mathsf{rec}}(E)$'
    xlabel = 'Energy (MeV)'
    ylabel = 'Particle emission yield per nt'
    #
    axes.set_title(title, y=position_title_y)
    axes.set_xlabel(xlabel, labelpad=labelpad_x)
    axes.set_ylabel(ylabel, labelpad=labelpad_y)
    axes.set_xlim([1.0E-3, Energy_max])
    axes.set_ylim([Yield_min, Yield_max])
    axes.set_xscale(x_scale)
    axes.set_yscale(y_scale)
    axes.legend(loc='best',  ncol=1, framealpha=0.95)
    axes.tick_params(axis='both', which='major')
    axes.grid(visible=True, which='major', axis='x',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='x',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.grid(visible=True, which='major', axis='y',
              color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
    axes.grid(visible=True, which='minor', axis='y',
              color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
    axes.minorticks_off()
    # fig.savefig(path_out_abs + nuclide + '_DER.svg', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='svg', transparent=True, pad_inches=0.00)
    # fig.savefig(path_out_abs + nuclide + '_DER.pdf', dpi=None, facecolor='white', edgecolor='black',
    #            orientation=None, format='pdf', transparent=True, pad_inches=0.00)
    # plt.show()
    plt.clf()
    plt.close()


# Plot ALL line-Spectra as stem-plots and ALL continuous spectra as line-plots on one axes:
fig, axes = plt.subplots(nrows=1, ncols=1, sharex=True, sharey=True)
fig.patch.set_linewidth(3)
fig.patch.set_edgecolor('xkcd:black')
#
if Yield_DEP != 0.0:
    x = eng_DEP
    y = Q_DEP
    (markers, stemlines, baseline) = plt.stem(
        x, y, label=r'$Q_{\gamma}(E)$', bottom=Yield_min)
    plt.setp(markers, marker='D', markersize=marker_size(x), markerfacecolor=color_markerface_DEP,
             markeredgecolor=color_markeredge_DEP, markeredgewidth=2, alpha=1.00)
    plt.setp(stemlines, linestyle='solid', color=color_stemline_DEP,
             linewidth=LineW_stem(x), alpha=1.00)
    plt.setp(baseline, linestyle='solid', color=color_baseline,
             linewidth=lw_baseline, visible=True, alpha=1.00)
#
if Yield_DEE != 0.0:
    x = eng_DEE
    y = Q_DEE
    (markers, stemlines, baseline) = plt.stem(
        x, y, label=r'$Q_{e}(E)$', bottom=Yield_min)
    plt.setp(markers, marker='o', markersize=marker_size(x), markerfacecolor=color_markerface_DEE,
             markeredgecolor=color_markeredge_DEE, markeredgewidth=2, alpha=1.0)
    plt.setp(stemlines, linestyle='solid', color=color_stemline_DEE,
             linewidth=LineW_stem(x), alpha=1.00)
    plt.setp(baseline, linestyle='solid', color=color_baseline,
             linewidth=lw_baseline, visible=True, alpha=1.00)
#
if Yield_DEA != 0.0:
    x = eng_DEA
    y = Q_DEA
    (markers, stemlines, baseline) = plt.stem(
        x, y, label=r'$Q_{\alpha}(E)$', bottom=Yield_min)
    plt.setp(markers, marker='o', markersize=marker_size(x), markerfacecolor=color_markerface_DEA,
             markeredgecolor=color_markeredge_DEA, markeredgewidth=2, alpha=1.00)
    plt.setp(stemlines, linestyle='solid', color=color_stemline_DEA,
             linewidth=LineW_stem(x), alpha=1.00)
    plt.setp(baseline, linestyle='solid', color=color_baseline,
             linewidth=lw_baseline, visible=True, alpha=1.00)
#
if Yield_DER != 0.0:
    x = eng_DER
    y = Q_DER
    (markers, stemlines, baseline) = plt.stem(
        x, y, label=r'$Q_{\mathsf{rec}}(E)$', bottom=Yield_min)
    plt.setp(markers, marker='h', markersize=marker_size(x), markerfacecolor=color_markerface_DER,
             markeredgecolor=color_markeredge_DER, markeredgewidth=2, alpha=1.00)
    plt.setp(stemlines, linestyle='solid', color=color_stemline_DER,
             linewidth=LineW_stem(x), alpha=1.00)
    plt.setp(baseline, linestyle='solid', color=color_baseline,
             linewidth=lw_baseline, visible=True, alpha=1.00)
#
if Yield_BEN != 0.0:
    x = eng_BEN
    y = Q_BEN
    axes.plot(x, y, label=r'$Q_{\beta^{-}}(E)$', color=color_BEN, linestyle='solid',
               marker=None, alpha=1.0)
#
if Yield_BEP != 0.0:
    x = eng_BEP
    y = Q_BEP
    axes.plot(x, y, label=r'$Q_{\beta^{+}}(E)$', color=color_BEP, linestyle='solid',
               marker=None, alpha=1.0)
#
if Yield_SFN != 0.0:
    x = eng_SFN
    y = Q_SFN
    axes.plot(x, y, label=r'$Q_\mathrm{{SFN}}(E)$', color=color_SFN, linestyle='solid',
               marker=None, alpha=1.00)
#
title  = LaTeX_RadioNuclide + r'  All-particle emission spectra $Q(E)$'
xlabel = 'Energy (MeV)'
ylabel_left = r'Discrete: Particle yield $\left(\mathsf{nt^{-1}} \right)$'
ylabel_right = r'Continuous: Particle yield $\left (\mathsf{nt^{-1} \; {MeV}^{-1}} \right)$'
axes_right = axes.twinx()
#
axes.set_title(title, y=position_title_y)
axes.set_xlabel(xlabel, labelpad=labelpad_x)
axes.set_ylabel(ylabel_left, labelpad=labelpad_y)
axes_right.set_ylabel(ylabel_right, labelpad=labelpad_y)
axes.set_xlim([1.0E-5, Energy_max])
axes.set_ylim([Yield_min, Yield_max])
axes.set_xscale(x_scale)  # 'log'
axes.set_yscale(y_scale)
axes_right.set_yscale(y_scale)
axes_right.set_ylim(axes.get_ylim())
axes.legend(loc='best',  ncol=1, framealpha=0.50)
axes.tick_params(axis='both', which='major')
axes_right.tick_params(axis='y', which='major')
axes.grid(visible=True, which='major', axis='x',
          color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
axes.grid(visible=True, which='minor', axis='x',
          color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
axes.grid(visible=True, which='major', axis='y',
          color=color_grid_major, alpha=1.0, linestyle='dashed', linewidth=lw_grid_major)
axes.grid(visible=True, which='minor', axis='y',
          color=color_grid_minor, alpha=1.0, linestyle='dotted', linewidth=lw_grid_minor)
axes.minorticks_off()
axes_right.minorticks_off()
fig.savefig(path_out_abs + nuclide + '_ALL.svg', dpi=None, facecolor='white', edgecolor='black',
            orientation=None, format='svg', transparent=True, pad_inches=0.00)
fig.savefig(path_out_abs + nuclide + '_ALL.pdf', dpi=None, facecolor='white', edgecolor='black',
            orientation=None, format='pdf', transparent=True, pad_inches=0.00)
plt.show()
plt.clf()
plt.close()

