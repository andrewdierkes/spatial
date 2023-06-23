#!/usr/bin/env python
# coding: utf-8

# <div class="alert alert-block alert-success">
# <b> Test results and statistical analysis of all cards sent by the screen printing team. This includes MESH 20, 31, 43, INK types A, E & I. With 3 bags of each Ink type (A/E) and two bags of I for each mesh type. <b> <\div>

# In[37]:


pip install pingouin


# In[38]:


pip install scikit_posthocs


# In[3]:


import re
import numpy as np
import pathlib as path
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy import stats
from datetime import date
import pingouin as pg
from itertools import combinations


# In[33]:


_path = path.Path.cwd()

assay_name_dataset = []

dataset_dict = {}
qos_dataset_dict = {}
spatial_list = []


for _filepath in _path.iterdir():
    if _filepath.suffix != r'.tsv':
        continue
    elif _filepath.suffix == r'.tsv':
        
        with open(_filepath, 'r') as file:
            
            #individual assay data across 4 channels
            assay_dict = {}
            wl_index = []
            
            data = file.read()
            
            assay_name_finder = re.compile(r'@DES:.*')
            pixelData_finder = re.compile(r'pixelData.*]')
            numberPixel_finder = re.compile(r'numberPixels=[0-9][0-9]')           
            beginPosition_finder = re.compile(r'beginPosition.*[0-9]')
            numberOfSteps_finder = re.compile(r'numberOfSteps.*[0-9]')
            stepSize_finder = re.compile(r'stepSize.*[0-9]')
            spd_finder = re.compile(r'COMMAND.*Sent:.*INS.*SPD.*')
            windowStartWavelength_finder = re.compile(r'windowStartWavelength.*[0-9] ')
            windowEndWavelength_finder = re.compile(r'Sent INS SSCPARAM.*[0-9]')
            qos_finder = re.compile(r'(?<=QOS Optical Signal: )[0-9]+.[0-9]+')          
            
            assay_name = re.findall(assay_name_finder, data)
            pixelData = re.findall(pixelData_finder, data)
            numberPixel = re.findall(numberPixel_finder, data)
            beginPosition = re.findall(beginPosition_finder, data)
            numberOfSteps = re.findall(numberOfSteps_finder, data)
            stepSize = re.findall(stepSize_finder, data)
            windowStartWavelength = re.findall(windowStartWavelength_finder, data)
            windowEndWavelength = re.findall(windowEndWavelength_finder, data)
            spd = re.findall(spd_finder, data)
            qos_value = qos_finder.findall(data)
        
            
            #find specific names of assay
            #print(file, assay_name)
            
            
            #number of scans you hope to take
            #scrub & transform data into float
            #if len(pixelData) == 80:    
            pixel_data_remove = re.compile(r'pixelData:.\[')
            pixel_data_removal = re.compile(r'\]')

            number_pixel_removal = re.compile(r'numberPixels=')
            begin_position_removal = re.compile(r'beginPosition.*=..')
            number_of_steps_removal = re.compile(r'numberOfSteps.*=..')
            step_size_removal = re.compile(r'stepSize.*=..')
            window_start_wl_removal = re.compile(r'windowStartWavelength=')
            spd_removal = re.compile('COMMAND\tSent: INS ')
            assay_name_remove = re.compile(r'@DES:.')


            pixel_data_semi = [pixel_data_remove.sub('', string) for string in pixelData]
            pixel_data = [var.split(', ') for var in [pixel_data_removal.sub('', string) for string in pixel_data_semi]]

            pixel_dataset = []

            for scan in pixel_data:
                scan_float = []
                for var in scan:
                    try:
                        scan_float.append(float(var))
                    except ValueError as e:
                            f'Cannot convert {var} to float'
                pixel_dataset.append(scan_float)
            
            

            #converting into floats and strings for parsing 
            number_pixel_list = [int(var) for var in [number_pixel_removal.sub('', string) for string in numberPixel]]
            #print(number_pixel_list, file)
            begin_position = [int(var) for var in [begin_position_removal.sub('', string) for string in beginPosition]]
            number_of_steps_list = [int(var) for var in [number_of_steps_removal.sub('', string) for string in numberOfSteps]]
            step_size_list = [int(var) for var in [step_size_removal.sub('', string) for string in stepSize]]
            window_start_wl_list = [float(var) for var in [window_start_wl_removal.sub('', string) for string in windowStartWavelength]]
            SPD = [spd_removal.sub('', string) for string in spd]
            window_end_wl_string = ''.join(windowEndWavelength)
            #ending WL edited as changing IT will effect the result (Old code: int(window_end_wl_string[40:43])
            window_end_wl = 650
            assay_list = [assay_name_remove.sub('', string) for string in assay_name]
                
            number_pixel = int(number_pixel_list[0])
            
            
            
            number_of_steps = int(number_of_steps_list[0])
            step_size = int(step_size_list[0])
            window_start_wl = int(window_start_wl_list[0])
            assay = assay_list[0]
            
            assay_name_dataset.append(assay)

            #using generator; divide scans into groups of 20 for each channel
            def chunks(lst, chunk):
                for i in range(0, len(lst), chunk):
                    yield lst[i:i+chunk]

            #unpack with generator into a list with *
            pixel_channels = [*chunks(pixel_dataset, 20)]
            pixel_ch0 = pixel_channels[0]
            pixel_ch1 = pixel_channels[1]
            pixel_ch2 = pixel_channels[2]
            pixel_ch3 = pixel_channels[3]
            

            #create spatial read values from data within assay
            def channel_chunk(ch_pos, step_size, num_steps):
                '''iterate thru spatial coords using ch_pos we start with, step size and num_steps from assay.. range(start,stop,interval)'''
                spatial_lister = []
                for i in range(ch_pos, ch_pos+(step_size*num_steps), step_size):
                                spatial_lister.append(i)
                return spatial_lister

            #spatial locations of assay used for lambda plotter...
            spatial_list.append(channel_chunk(int(begin_position[0]), step_size, number_of_steps))
            spatial_list.append(channel_chunk(int(begin_position[1]), step_size, number_of_steps))
            spatial_list.append(channel_chunk(int(begin_position[2]), step_size, number_of_steps))
            spatial_list.append(channel_chunk(int(begin_position[3]), step_size, number_of_steps))

            #create index using starting and ending wl and number of pixels in each scan array; adding outside of loop
            wl_step = ((window_end_wl-window_start_wl)/number_pixel)
            start_wl = (window_start_wl-wl_step)
        
            #ending WL edited as changing IT will effect the result (Old code=int(window_end_wl_string[40:43])-wl_step)
            end_wl = 650-wl_step

            offset = start_wl 

            while offset < end_wl:

                i = offset + wl_step
                wl_index.append(i)
                offset += wl_step

                #add spectral and spatial data to outside of log reading script
            assay_dict[f'{assay} ch0'] = dict(zip(spatial_list[0], pixel_ch0))
            assay_dict[f'{assay} ch1'] = dict(zip(spatial_list[1], pixel_ch1))
            assay_dict[f'{assay} ch2'] = dict(zip(spatial_list[2], pixel_ch2))
            assay_dict[f'{assay} ch3'] = dict(zip(spatial_list[3], pixel_ch3))

            dataset_dict[f'{assay}'] = assay_dict
            
                        
            #qos data
            qos_dict = {}
            qos_dict['assay_name'] = assay
            qos_dict['channels'] = (1400,2200,3000,3800)
            qos_dict['qos_value'] = [float(var) for var in qos_value]
            
            #find assays which are not correct
            for var in qos_dict['qos_value']:
                if var < 10000:
                    print(file,var)
                    
            qos_dataset_dict[f'{assay}'] = qos_dict

            pd.set_option('display.max_columns', None)
            pd.set_option('display.max_rows', None)      
          

            #else:
             #   print(f'number of scans if off from 80, check file {file}')


# In[34]:


def delist(args):
    '''delists a list of lists into 1 list'''
    delist = [var for small_list in args for var in small_list]
    return(delist) 


# In[6]:


#lambda plotter will be used in an iterator below for all assays
data_dict = {}

def lambda_plotter(df, wl_index, channel_sp, graph_channel):
    '''df = spatial dataframe, wl_index = 520-650, 
    channel_sp = spatial channel data, graph_channel = assay_name'''
    
    def delist(args):
        '''delists a list of lists into 1 list'''
        delist = [var for small_list in args for var in small_list]
        return(delist) 

    ch_max = [df.max(axis=1)]
    ch_spatial = channel_sp
    
    #iterate thru max list 
    for var in ch_max:
        lambda_max = (max(var))

    lmax_index = []

    #iterate thru df columns looking for the lamda max value
    for col in df.columns:
        index_lookup = df[df[col] == lambda_max].index.tolist()
        try:
        #if a variable exists... appended
            lmax_index.append(index_lookup[0])
        except:
            pass

    var = float(lmax_index[0])
    
    lambda_max_row = wl_index.index(var)
    
    #avg function
    def average(data, length):
        avg = sum(data)/len(data)
        return avg
    
    #create an average list similar to plateau if it is within 10% of the RFU
    RFU_max_list = []
    for z in df.iloc[lambda_max_row]:
        if z >= lambda_max*.92:
            RFU_max_list.append(z)
    
    num_spatial_pos = len(RFU_max_list)
    
    average = average(RFU_max_list, num_spatial_pos)

    
    #take the spatial data across the max wavelength (df.iloc) and assign variables as arrays for find_peaks
    x = np.array(channel_sp)
    y = np.array(df.iloc[lambda_max_row])
    
    #CURRENTLY NOT USING
    #scipy peak finder; use 1d np.array and height = min value for peak 
    peaks = find_peaks(y, height = 2000, distance=10)
    peak_pos = x[peaks[0]]
    height = peaks[1]['peak_heights']
    
    
    
    #plot spatial vs RFU at max wl & maxima via scipy
#     fig, ax = plt.subplots()
#     plt.scatter(ch_spatial, y)
#     plt.scatter(peak_pos, height, color = 'r', s = 30, marker = 'D', label = 'maxima')
    

#     ax.set(xlabel='Spatial Position', ylabel='RFU', title= f'Formulation {graph_channel} Spatial RFU & Lambda Max at {round(var, 3)}')
#     ax.legend()
#     ax.grid()
    
    #create columns from assay_name for statistical sorting
    mesh_finder = re.compile(r'(?<=M)[0-9]+')
    ink_finder = re.compile(r'(?<=I)[A-Z]+')
    split_finder = re.compile(r'(?<=S)[0-9]+')
    card_finder = re.compile(r'(?<=C)[0-9]+')

    mesh = mesh_finder.findall(graph_channel)
    ink = ink_finder.findall(graph_channel)
    split = split_finder.findall(graph_channel)
    card = card_finder.findall(graph_channel)
    mesh_ink = '_'.join([mesh[0],ink[0]])
    
    #data dictionary for each assay
    data_dict['assay_name'] = graph_channel
    data_dict['mesh_ink'] = mesh_ink
    data_dict['mesh'] = int(mesh[0])
    data_dict['ink'] = ink[0]
    data_dict['split'] = split[0]
    data_dict['card'] = card[0]
    data_dict['max_rfu'] = lambda_max #height
    data_dict['average_8'] = average 
    data_dict['num_spatial_pos_avg'] = num_spatial_pos 
    data_dict['lambda_max'] = var
    #data_dict['actual_spatial_pos_max_rfu'] = peak_pos
    
    #return f'Max RFU output {lambda_max}, seen at wavelength: {var}'


# In[7]:


def heat_map(dataframe):
    hm = sns.heatmap(dataframe, cmap='YlOrRd')


# In[8]:


#generate a list of the assay names pairing number of repeats with number of channels
assay_single = []
for assay in dataset_dict.keys():
    assay_single.append(np.repeat(assay, 4))
    
assay_list = delist(assay_single)

#for lambda plotter title
channel = [1400,2200,3000,3800]
channel_assay = channel*(len(dataset_dict.keys()))


#create an iterator to go thru assay_list and pair the assay_name with the four corresponding dataframes
assay_chunk = 1
offset = 0


summary_list = []

#iterate thru dataset_dict getting into specific channels to create dataframes & link assay_name with DF
for assay in dataset_dict.values():
    for channel in assay.values():
        
        df = pd.DataFrame(channel, index=wl_index)
        df.index.name = assay_list[offset]
        
        
        #list_of_data.append(df)
        
        #display(df)
        lambda_plotter(df, wl_index, spatial_list[offset], assay_list[offset])
        df_summary = pd.DataFrame(data_dict, index=[0])
        summary_list.append(df_summary)
        
        offset += assay_chunk
        
        #display(heat_map(df))
        
        #display(hm = sns.heatmap(df, cmap='YlOrRd'))
        


# In[12]:


#extract data on QOS values
qos_summary = []
for key, value in qos_dataset_dict.items():
    df_qos = pd.DataFrame(value)
    qos_summary.append(df_qos)


# <div class="alert alert-block alert-warning">
# <b>First, we'll take a look at the spatial data extracted from the meter: <b> 
# </div>

# In[13]:


def chunk_cv(df_col, chunk=4):
    '''This function will iterate over a df_col and return the average and stdev, using chunk_iterator.. so if you'd like the average & stdev values occuring every 3 rows... use 3 as your chunk iterator'''

    offset_mean = 0
    offset_stdev = 0
    
    number_list = [var for var in range(len(df_col))]
    
    dataset_average = []
    dataset_stdev = []
    
    
    while offset_mean < len(number_list):
        i_mean = number_list[offset_mean:chunk+offset_mean]
        average = df_col.iloc[i_mean].mean(axis=0)
        
        dataset_average.append(average)
        #dataset_array.append(_array)
        
        offset_mean += chunk
    
    while offset_stdev < len(number_list):
        i_stdev = number_list[offset_stdev:chunk+offset_stdev]
        stdev = df_col.iloc[i_stdev].std(ddof=1)
        
        dataset_stdev.append(stdev)
        
        offset_stdev += chunk
    
    return dataset_average, dataset_stdev


# In[14]:


def unique(df_col, chunk=4):
    '''This function will iterate over a df_col and return only unique values using chunk_iterator.. so if you'd like the unique value of something occuring every 3 rows... use 3 as your chunk iterator'''

    offset = 0
    
    number_list = [var for var in range(len(df_col))]
    
    dataset_array = []
    dataset_list = []
    
    while offset < len(number_list):
        i = number_list[offset:chunk+offset]
        _array = df_col.iloc[i].unique()
        
        dataset_array.append(_array)
        
        offset += chunk
    
    for var in dataset_array:
        dataset_list.append(var.tolist())
    
    unique = delist(dataset_list)
    
    return unique


# In[15]:


#find unique names from Assay Name row, slice each name to get rid of the rep number and then pass through unique function
#name_len = length of names we want to slice

def slice_name(df, row, name_len, end_offset):
    '''Function SLICE_NAME slices strings (usually assay name) that are similar in nature (df row). It works to 
    remove the ends (end_offset) which are different (usually rep number) to allow the UNIQUE function to work
    df = dataframe to use
    row = the index of the column you want to slice
    name_len = the length of each assay name, they all should be the same
    end_offset = what position you want to end with'''
    
    assay__name = [var for var in df.iloc[:,row]]

    unique_name = []


    for var in assay__name:
        if len(var) == name_len:
            unique_name.append(var[0:end_offset])

        else:
            print(file)
            print(f'labeling error for assay:{var}, length {len(var)}')
            for count, var2 in enumerate(var):
                print(count, var2)
    
    return pd.DataFrame(unique_name)


# <div class="alert alert-block alert-warning">
# <b>Below is a dataframe & associated compilation of boxplots holding average spatial RFU reads and their associated intrastrip cv <b> 
# </div>

# In[78]:


#raw spatial data


df_dataset_sum = pd.concat(summary_list, axis=0).sort_index().reset_index(drop=True)
df_dataset_sum.sort_values(by=['assay_name'], inplace=True)

#add in target spatial_positions
ch = (1400,2200,3000,3800)
repeat_ch =(ch*int((len(df_dataset_sum)/4)))

df_dataset_sum.insert(10, 'target_spatial_pos_max_rfu', repeat_ch)

print(len(df_dataset_sum))
#display(df_dataset_sum)
#print(len(df_dataset_sum))
df_dataset_sum.boxplot(column='max_rfu', by='mesh_ink', figsize= (8,6), ylabel='RFU')
#df_dataset_sum.boxplot(column='max_rfu', by='ink', figsize= (8,6))


# In[49]:


#find assays with more than 4 channels recorded
assay_list = []
for var in df_dataset_sum.iloc[:,0]:
    assay_list.append(var)

for var in df_dataset_sum.iloc[:,0]:
    find_assay = re.compile(f'{var}')
    num_ch = find_assay.findall(''.join(assay_list))
    if len(num_ch) != 4:
        print(num_ch)



# <div class="alert alert-block alert-warning">
# <b> Next, well take a look at spatial data looking at intrastrip values such as AVG RFU, Varience & CV  (average of four channels from each strip) & box plots comparing mesh, ink and mesh_ink categories <b>
# </div>

# In[24]:


#intrastrip based on spatial data
s_avg, s_stdev = chunk_cv(df_dataset_sum.iloc[:,6],4)
s_name = unique(df_dataset_sum.iloc[:,0], 4) 

spat_strip_dict = {}
spat_strip_dict['assay'] = s_name
spat_strip_dict['average_rfu'] = s_avg
spat_strip_dict['stdev'] = s_stdev

df_spatial_strip = pd.DataFrame(spat_strip_dict)

#put metadata back into dataframe
md_mesh = []
md_ink = []
md_split = []
md_card = []
md_meshink = []

for var in df_spatial_strip.iloc[:,0]:
    
    mesh_finder = re.compile(r'(?<=M)[0-9]+')
    ink_finder = re.compile(r'(?<=I)[A-Z]+')
    split_finder = re.compile(r'(?<=S)[0-9]+')
    card_finder = re.compile(r'(?<=C)[0-9]+')

    mesh = mesh_finder.findall(var)
    ink = ink_finder.findall(var)
    split = split_finder.findall(var)
    card = card_finder.findall(var)
    mesh_ink = '_'.join([mesh[0],ink[0]])
    
    
    md_mesh.append(int(''.join(mesh)))
    md_ink.append(''.join(ink))
    md_split.append(int(''.join(split)))
    md_card.append(int(''.join(card)))
    md_meshink.append(''.join(mesh_ink))



#insert important metadata back into df
df_spatial_strip.insert(1, 'card', md_card)
df_spatial_strip.insert(1, 'split', md_split)
df_spatial_strip.insert(1, 'ink', md_ink)
df_spatial_strip.insert(1, 'mesh', md_mesh)
df_spatial_strip.insert(1, 'mesh_ink', md_meshink)
df_spatial_strip.insert(8, 'intrastrip_cv', round(((df_spatial_strip.iloc[:,7]/df_spatial_strip.iloc[:,6])*100),4))

print(len(df_spatial_strip))
display(df_spatial_strip.head())

#print(f'Average intrastrip CV based on spatial data is {df_spatial_strip.iloc[:,3].mean()}')


display(df_spatial_strip.boxplot(column='intrastrip_cv', by='mesh_ink', figsize= (8,6), ylabel='intrastrip_cv'))
display(df_spatial_strip.boxplot(column='intrastrip_cv', by='ink', figsize= (8,6), ylabel='intrastrip_cv'))
display(df_spatial_strip.boxplot(column='intrastrip_cv', by='mesh', figsize=(8,6), ylabel='intrastrip_cv'))


# In[133]:


def regex_cv(df,groupby,feature):
    '''regex_cv applies the inputted patterns on df_feature.
    df = dataframe, groupby = col # for where to iloc and search for the regex,
    feature = col # which you want to aggregate mean, stdev & cv'''
    
    regex = df.iloc[:,groupby].unique()
    
    groupby_name = df.columns[groupby]
    
    #sort dataframe to align with regex iterator below(min,max... basically need groupings in one space)
    df.sort_values(by=[groupby_name], inplace=True)
    df.reset_index(drop=True, inplace=True)

    display(df.head(), df.tail())
    
    print(f'grouped by {groupby_name}')
    #find indexes of regex patterns
    idx_list = []
    for var in regex:
        idx = df[df.iloc[:,groupby] == var].index.to_list()
        idx_list.append(idx)

    #iterate thru df using idx_list
    cv_list = []
    mean_list = []
    stdev_list = []
    
    for var in idx_list:
        #print(var)
        if len(var) == 0:
            pass
        else:
            _min = min(var)
            _max = max(var)
            #print(_min,_max)
            chunk = len(var)
            offset = 0

            mean = df.iloc[_min:_max, feature].mean()
            stdev = df.iloc[_min:_max, feature].std()
            cv = (stdev/mean)*100
    
            
            mean_list.append(mean)
            stdev_list.append(stdev)
            cv_list.append(cv)
            
    final = list(zip(cv_list,mean_list,stdev_list))
    df_final = pd.DataFrame(final, index=[regex], columns=['cv','mean','stdev'])
    display(df_final)
    
    return df_final


# In[118]:


def normality(df,groupby,feature):
    '''df=dataframe name, groupby=different categories were comparing(IV),
    feature=DV,'''
    
    import pylab
    
    regex = df.iloc[:,groupby].unique()
    
    groupby_name = df.columns[groupby]
    
    #sort dataframe to align with regex iterator below(min,max... basically need groupings in one space)
    df.sort_values(by=[groupby_name], inplace=True)
    df.reset_index(drop=True, inplace=True)

    display(df.head(),df.tail())
    
    print(f'grouped by {groupby_name}')
    #find indexes of regex patterns
    idx_list = []
    for var in regex:
        idx = df[df.iloc[:,groupby] == var].index.to_list()
        idx_list.append(idx)

    #iterate thru df using idx_list
    cv_list = []
    i =0
    for var in idx_list:
        if len(var) == 0:
            pass
        else:
            _min = min(var)
            _max = max(var)
            chunk = len(var)
            offset = 0
            print(_min,_max,type(_min))
                    #KDE & QQ plots
            plt.figure(figsize=(10,5))
            plt.subplot(1,2,1, title=f"Grouping {''.join(str(regex[i]))}")
            sns.kdeplot(df.iloc[_min:_max,feature])
            plt.subplot(1,2,2)
            stats.probplot(df.iloc[_min:_max,feature],plot=pylab)
            plt.show()
            i += 1


# <div class="alert alert-block alert-warning">
# <b> Statistical analysis of different mesh types. Below we have KDE plots and Q-Q plot (plot 2 sets of quantiles against each other. The Kernel Density Estimator reveals the distribution of samples in a simple curve. 
#     A quantile is a division of the data sorted in ascending order. We compare a quantile of our data against a theretical quantile with normal distribution. Looking for a slope of 1.
#     We also run a shapiro wilks hypothesis test for normality, rejecting the null (groups distribution is guassian) as p < 0.05 <b>
# </div>

# In[119]:


#comparing by mesh

display(df_spatial_strip.boxplot(column='intrastrip_cv', by='mesh', figsize= (7,5), ylabel='intrastrip_cv'))


#kde plot & qq plot
normality(df_spatial_strip,2,8) 

#check distribution 
print('shapiro-wilks test for gaussian distribution; p-val < 0.05; reject null that groups are gaussian.')
pg.normality(df_spatial_strip, dv='intrastrip_cv', group='mesh', method='shapiro')



# <div class="alert alert-block alert-warning">
# <b> Below we group the intrastrip CVs and find the average CV per mesh grouping; as of now we like mesh 31" <b>
# </div>

# In[120]:


mesh_cvavg, mesh_cvstdev = chunk_cv(df_spatial_strip.iloc[:,8],160)
df_mesh_avg = pd.DataFrame(mesh_cvavg, columns=['mesh_avg_cv'],index=[df_spatial_strip.iloc[:,2].unique()])
display(df_mesh_avg)


# <div class="alert alert-block alert-warning">
# <b> shapiro-wilks hypothesis tests are all false. This means we cannot do an ANOVA as it assumes normal distribution and varience among groups. We perform a Kruskal Wallis non-parametric test to determine if there are significant differences in the medians of each group <b>
# </div>

# In[121]:


#kruskal-wallis
def kruskal_array(df, groupby, feature):
    '''df = dataframe, groupby=iv,
    feature = dv'''
    regex = df.iloc[:,groupby].unique()
    
    groupby_name = df.columns[groupby]
    
    #sort dataframe to align with regex iterator below(min,max... basically need groupings in one space)
    df.sort_values(by=[groupby_name], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    #display(df)
    display(df.head(), df.tail())
    
    print(f'grouped by {groupby_name}')
    
    #find indexes of regex patterns
    idx_list = []
    for var in regex:
        idx = df[df.iloc[:,groupby] == var].index.to_list()
        idx_list.append(idx)
    
    for var in idx_list:
        print(f'length of group = {len(var)}')
        
    #iterate thru df using idx_list
    kw_list = []
    kw_dict = {}
    for var in idx_list:
        if len(var) == 0:
            pass
        else:
            _min = min(var)
            _max = max(var)
            chunk = len(var)
            offset = 0
            print(_min,_max,chunk)
            kw_list.append(np.array(df.iloc[_min:_max, feature]))
            kw_dict[f'{df.iloc[offset, groupby]}'] = np.array(df.iloc[_min:_max, feature])
          
    return kw_list
        


# <div class="alert alert-block alert-warning">
# <b> Next we take a look at the KW hypothesis test seeing if there is a difference in medians for the mesh categories <b>
# </div>

# In[122]:


kw_mesh = kruskal_array(df_spatial_strip, 2, 8)

result = stats.kruskal(kw_mesh[0],kw_mesh[1], kw_mesh[2])
print(result)

if result[1] < 0.05:
    print(f'our p-val of {result[1]} is less than our significance value of 0.05. This means we can reject our null & know that the medians of k groups differ')
if result[1] > 0.05:
     print(f'our p-val of {result[1]} is greater than our significance value of 0.05. This means we cant our null & know that the medians of k groups come from similar populations')


# <div class="alert alert-block alert-warning">
# <b> We get a H value of 19.29 and a pval < 0.05. Because we have a low chance of the medians being from the same population distribution; we reject the null hypothesis, we know there are groups which differ and we'd like to see which groups differ from others so we run the dunn post hoc with a bonferroni p-value adjustment (to reduce family wise error) <b>
# </div>

# In[123]:


#POST HOC for non parametric test: Drunns with a bonferroni correction to reduce family wise error
import scikit_posthocs as sp

mesh_dunn = sp.posthoc_dunn(kw_mesh, p_adjust = 'bonferroni')
display(mesh_dunn)


#if p-val < 0.05 -> reject null and disprove likeness... difference in groups median
print('significant difference between groups 1 & 2 as well as 3 & 2. Proving 2 is the odd one out')


# <div class="alert alert-block alert-warning">
# <b> Below we have the data sorted for strip to strip variation in each mesh-ink type (number of cards = 3 A/E & 3 for I). The top two are df head/tail & the bottom is the CV, mean & stdev for each mesh_ink type (n=60 a/e, n=40 i)<b>
# </div>

# In[124]:


regex_cv(df_spatial_strip, 1,6)


# <div class="alert alert-block alert-success">
# <b> Statistical analysis of different ink types <b>
# </div>

# In[125]:


df_spatial_strip.boxplot(column='intrastrip_cv',by='ink', figsize=(8,6), ylabel='cv')


# <div class="alert alert-block alert-warning">
# <b>Looking at the strip to strip averages, variences, and cvs. Below we are comparing all cards across different meshes of the same ink type</b>

# In[126]:


print(normality.__doc__)
df_spatial_strip.head()
normality(df_spatial_strip, 3, 8)

pg.normality(df_spatial_strip, dv='intrastrip_cv', group='ink', method='shapiro')


# <div class="alert alert-block alert-warning">
# <b> Looking at the KDE and q-q plot once again, with a new perspective (INK) and running a shapiro-wilks hypothesis test; we see none of the inks follow a gaussian distribution, so we choose a Kruskal Wallis analysis" <b>
# </div>

# In[127]:


#non parametric kruskal 
kruskal_array.__doc__
kw_ink = kruskal_array(df_spatial_strip, 3, 8)


result = stats.kruskal(kw_ink[0],kw_ink[1], kw_ink[2])
print(result)

if result[1] < 0.05:
    print(f'our p-val of {result[1]} is less than our significance value of 0.05. This means we can reject our null & know that the medians of k groups differ')
if result[1] > 0.05:
     print(f'our p-val of {result[1]} is greater than our significance value of 0.05. This means we cant our null & know that the medians of k groups come from similar populations')


# <div class="alert alert-block alert-warning">
# <b> With a H value of 5.35 & a p-val of 0.06, we cant reject our null; meaning that the inks A, E, and I come from similar population distributions with medians that are similar; we will not run a post-hoc dunn's anaylsis
#     <b>
# </div>

# <div class="alert alert-block alert-warning">
# <b> Here we take a look at the differences in channel RFU across the same mesh_ink type
#     <b>
# </div>

# In[86]:


#Compare Channels of similar position using spatial data
#channel_comparison

df_dataset_sum.sort_values(by=['assay_name'])
#display(df_dataset_sum)


ch_1400 = [var for var in range(0,len(df_dataset_sum), 4)]
ch_2200 = [var for var in range(1, len(df_dataset_sum),4)]
ch_3000 = [var for var in range(2, len(df_dataset_sum),4)]
ch_3800 = [var for var in range(3, len(df_dataset_sum),4)]



ch_1400_df = []
for var in ch_1400:
    _1400= df_dataset_sum.iloc[var]
    ch_1400_df.append(_1400)
df_ch1400 = pd.DataFrame(ch_1400_df)
#display(df_ch1400)

ch_2200_df = []
for var in ch_2200:
    _2200= df_dataset_sum.iloc[var]
    ch_2200_df.append(_2200)
df_ch2200 = pd.DataFrame(ch_2200_df)
#display(df_ch2200)

ch_3000_df = []
for var in ch_3000:
    _3000= df_dataset_sum.iloc[var]
    ch_3000_df.append(_3000)
df_ch3000 = pd.DataFrame(ch_3000_df)
#display(df_ch3000)

ch_3800_df = []
for var in ch_3800:
    _3800= df_dataset_sum.iloc[var]
    ch_3800_df.append(_3800)
df_ch3800 = pd.DataFrame(ch_3800_df)
#display(df_ch3800)


# In[174]:


#channel to channel cv based on spatial reads
print('Channel strip to strip CVs starting at 1400 and ending with 3800')
df_channel_dataset = pd.concat(
                            [df_ch1400.reset_index(drop=True), 
                             df_ch2200.reset_index(drop=True), 
                             df_ch3000.reset_index(drop=True),
                             df_ch3800.reset_index(drop=True)])
#display(df_channel_dataset)
                    


spatial_stos = {}

spatial_stos['cv_ch1400'] = regex_cv(df_ch1400, 1, 6)
spatial_stos['cv_ch2200'] = regex_cv(df_ch2200, 1, 6)
spatial_stos['cv_ch3000'] = regex_cv(df_ch3000, 1, 6)
spatial_stos['cv_ch3800'] = regex_cv(df_ch3800, 1, 6)


#df_spatial_stos = pd.DataFrame(spatial_stos, index=[df_ch1400.iloc[:,1].unique()])
#display(df_spatial_stos)
print(len(df_channel_dataset))

#print(f'Average strip to strip channel cvs for spatial data is {df_spatial_stos.iloc[0].mean()} cv')


# <div class="alert alert-block alert-warning"> <b> Next, we are evaluate the channel to channel variation seen amongst the different mesh_ink types. It is evident A across all mesh types has a lower ch-ch CV across all channels: </div>

# In[175]:


for key,value in spatial_stos.items():
    value.sort_values(by=['cv'],inplace=True)
    value.index.set_names(key, inplace=True)
    display(value)


# In[ ]:




