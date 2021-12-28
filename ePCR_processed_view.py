#!/usr/bin/env python
# coding: utf-8

# In[ ]:

from __future__ import division
import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from urllib.error import URLError
import math
import numpy as np
from pandas import datetime
import warnings
from io import StringIO
warnings.filterwarnings('ignore')
from typing import Callable
pd.options.display.float_format = '{:,.0f}'.format
from pathlib import Path



st.set_page_config(layout="wide")

st.title('ePCR viewer')

st.subheader("Upload already parsed data from ePCR - use below to view total period data while below individual arrays can be viewed. Uploading of data size greater than 200 Mb is not recommended")

st.button("Re-load")



uploaded_files = st.file_uploader("Choose a processed ePCR file ", accept_multiple_files=True)
for uploaded_file in uploaded_files:
    bytes_data = uploaded_file.read()
    bytes_data = str(bytes_data , "UTF-8")
   # st.write("filename:", uploaded_file.name)
    data = StringIO(bytes_data)
    #st.write(data) 
    print(type(data))
    comp = pd.read_csv(data, sep=",")
    conditions = [
        (comp['norm_N_Cov'] <= 4.0) & (comp['norm_RNaseP'] > 2.0),
        (comp['norm_N_Cov'] > 4.0) & (comp['norm_N_Cov'] <= 9.0) & (comp['norm_RNaseP'] >1.0),
        (comp['norm_N_Cov'] >= 9.0) & (comp['norm_RNaseP'] >=1.0),
        (comp['norm_N_Cov'] >= 9.0) & (comp['norm_RNaseP']<= 1.0),
        (comp['norm_N_Cov'] <= 4.0) & (comp['norm_RNaseP'] <=2.0),
        (comp['norm_N_Cov'] > 3.0) & (comp['norm_N_Cov'] <= 9.0) & (comp['norm_RNaseP'] <1.0)]

# create a list of the values we want to assign for each condition
    values = ['Negative_sample', 'PLOD', 'N_Cov_Positive_Sample', 'Control_N_Cov', 'No_Call','Background_PLOD']

     # create a new column and use np.select to assign values to it using our lists as arguments
    comp['Result'] = np.select(conditions, values)




#percentiles
def Q25(x):
    return x.quantile(0.25)

def Q50(x):
    return x.quantile(0.5)

def Q75(x):
    return x.quantile(0.75)


def ROXCV(df1):
    stats_ROX = df1.groupby(['Run_ID'])['ROX_RFU'].agg(['count', 'mean','std','min',Q25, Q50, Q75, 'max'])
    print('-'*30)
    CI95_hi_ROX = []
    CI95_lo_ROX = []
    CV_run_ROX = []
    for i in stats_ROX.index:
        c,m,s,t,u,q1,q2,v =(stats_ROX.loc[i])
        CI95_hi_ROX.append(m + 1.95*s/math.sqrt(c))
        CI95_lo_ROX.append(m - 1.95*s/math.sqrt(c))
        CV_run_ROX.append(s/m*100)
        
    stats_ROX['CI95% low ROX'] = CI95_lo_ROX
    stats_ROX['CI95% hi ROX'] = CI95_hi_ROX
    stats_ROX['ROX CV%'] = CV_run_ROX
    stats_ROX = stats_ROX.reset_index()
    return(stats_ROX)

stats_ROX = ROXCV(comp)



def fam_stat(df):
    stats_FAM = df.groupby(['Run_ID','Result'])['FAM_RFU'].agg(['count', 'mean','std','min',Q25, Q50, Q75, 'max'])
    print('stats_FAM')
    print('-'*30)
    
    CI95_hi_FAM = []
    CI95_lo_FAM = []
    CV_ruFAM_RFU = []
    for i in stats_FAM.index:
        c,m,s,t,u,q1,q2,v = stats_FAM.loc[i]
        CI95_hi_FAM.append(m + 1.95*s/math.sqrt(c))
        CI95_lo_FAM.append(m - 1.95*s/math.sqrt(c))
        CV_ruFAM_RFU.append((s/m*100))

    stats_FAM['ci95_lo_FAM'] = CI95_lo_FAM
    stats_FAM['ci95_hi_FAM'] = CI95_hi_FAM
    stats_FAM['CV%_FAM'] = CV_ruFAM_RFU
    stats_FAM = stats_FAM.reset_index().fillna('-')
    return(stats_FAM)


stats_FAM = fam_stat(comp)

def nfam_stat(df):
    stats_nFAM = df.groupby(['Run_ID','Result'])['norm_N_Cov'].agg(['count', 'mean','std','min',Q25, Q50, Q75, 'max'])
    print('stats_nFAM')
    print('-'*30)
    
    CI95_hi_nFAM = []
    CI95_lo_nFAM = []
    CV_nFAM_RFU = []
    for i in stats_nFAM.index:
        c,m,s,t,u,q1,q2,v = stats_nFAM.loc[i]
        CI95_hi_nFAM.append(m + 1.95*s/math.sqrt(c))
        CI95_lo_nFAM.append(m - 1.95*s/math.sqrt(c))
        CV_nFAM_RFU.append((s/m*100))

    stats_nFAM['ci95_lo_nFAM'] = CI95_lo_nFAM
    stats_nFAM['ci95_hi_nFAM'] = CI95_hi_nFAM
    stats_nFAM['CV%_nFAM'] = CV_nFAM_RFU
    stats_nFAM = stats_nFAM.reset_index().fillna('-')
    return(stats_nFAM)


stats_nFAM = nfam_stat(comp)

fig2b = px.scatter(comp, x= 'norm_RNaseP', y = 'norm_N_Cov',color = 'Result')
fig2b.update_traces(marker_size=3)

#fig2b.show()
#fig2b.write_html('comp_N3.html')


fig1bbnbb = px.scatter(comp, x= 'order', y = 'norm_RNaseP', color = 'Result')
fig1bbnbb.update_traces(marker_size=3)
fig1bbnbb.update_yaxes(range=[0, 6])
fig1bbnbb.add_trace(go.Scatter(
     y=[1.5, 1.5],
     x=[comp.order.min(), comp.order.max()],
     mode="lines+markers+text",
     name="Lower_1.5_RP Detected_Boundary",
     text=["1.5"],
     textposition="top center"))
fig1bbnbb.add_trace(go.Scatter(
     y=[2, 2],
     x=[comp.order.min(), comp.order.max()],
     mode="lines+markers+text",
     name="Lower_2_RP Detected_Boundary",
     text=["2"],
     textposition="top center"))
#fig1bbnbb.show()
#fig1bbnbb.write_html('comp_N3__monitor_normRNaseP.html')

figROX = px.scatter(comp, x= 'order', y = 'ROX_RFU', color = 'Result', title = 'Dispense Trace ROX')
figROX.update_yaxes(range=[1000, 6000], gridwidth = 0.0002, gridcolor ='grey')
figROX.update_traces(marker_size=3)

figROX.add_trace(go.Scatter(
    x=[comp.order.min(), comp.order.max()],
    y=[1600, 1600],
    mode="lines",
    name="1600  RFU Lower Cutoff Limit",
    text=["LCL"],
    #text=["ROX 1600 lower cutoff"],
    textposition="top center",
    line=dict(color="grey")
))
#figrp.show()


figrp = px.scatter(comp, x= 'ROX_RFU', y = 'FAM_RFU' ,color = 'Result')
figrp.update_traces(marker_size=3)
figrp.update_xaxes(range=[1000, 6000])
figrp.update_yaxes(range=[0, 50000])


figrp.add_trace(go.Scatter(
    x=[1600, 1600],
    y=[50000, -100],
    mode="lines",
    name="1600  RFU Lower Cutoff Limit",
    text=["LCL"],
    #text=["ROX 1600 lower cutoff"],
    textposition="top center",
    line=dict(color="grey")
))
#figrp.show()

def CFO_stat(df):
    stats_CFO = df.groupby(['Run_ID','Result'])['VIC_RFU'].agg(['count', 'mean','std','min',Q25, Q50, Q75, 'max'])
    print('stats_nFAM')
    print('-'*30)
    
    CI95_hi_CFO = []
    CI95_lo_CFO = []
    CV_CFO_RFU = []
    for i in stats_CFO.index:
        c,m,s,t,u,q1,q2,v = stats_CFO.loc[i]
        CI95_hi_CFO.append(m + 1.95*s/math.sqrt(c))
        CI95_lo_CFO.append(m - 1.95*s/math.sqrt(c))
        CV_CFO_RFU.append((s/m*100))

    stats_CFO['ci95_lo_CFO'] = CI95_lo_CFO
    stats_CFO['ci95_hi_CFO'] = CI95_hi_CFO
    stats_CFO['CV%_CFO'] = CV_CFO_RFU
    stats_CFO = stats_CFO.reset_index().fillna('-')
    return(stats_CFO)

stats_CFO = CFO_stat(comp)

def nCFO_stat(df):
    stats_nCFO = df.groupby(['Run_ID','Result'])['norm_RNaseP'].agg(['count', 'mean','std','min',Q25, Q50, Q75, 'max'])
    print('stats_nCFO')
    print('-'*30)
    
    CI95_hi_nCFO = []
    CI95_lo_nCFO = []
    CV_nCFO_RFU = []
    for i in stats_nCFO.index:
        c,m,s,t,u,q1,q2,v = stats_nCFO.loc[i]
        CI95_hi_nCFO.append(m + 1.95*s/math.sqrt(c))
        CI95_lo_nCFO.append(m - 1.95*s/math.sqrt(c))
        CV_nCFO_RFU.append((s/m*100))

    stats_nCFO['ci95_lo_CFO'] = CI95_lo_nCFO
    stats_nCFO['ci95_hi_CFO'] = CI95_hi_nCFO
    stats_nCFO['CV%_CFO'] = CV_nCFO_RFU
    stats_nCFO = stats_nCFO.reset_index().fillna('-')
    return(stats_nCFO)

stats_nCFO = nCFO_stat(comp)


figN1 = px.scatter(comp, x= 'order', y = 'norm_N_Cov' ,color = 'Result', title = 'N1 N2 Calls')

figN1.add_trace(go.Scatter(
    y=[10, 10],
    x=[comp.order.min(), comp.order.max()],
    mode="lines+markers+text",
    name="Lower_10_Positive_Boundary",
    text=["10"],
    textposition="top center",
    line=dict(color="red")
))

figN1.update_traces(marker_size=3)

figN1.add_trace(go.Scatter(
     y=[9, 9],
     x=[comp.order.min(), comp.order.max()],
     mode="lines+markers+text",
     name="Lower_9_Positive_Boundary",
     text=["9"],
     textposition="top center"))


figN1.update_xaxes(showgrid = True, gridwidth = 0.0002, gridcolor = 'grey')
figN1.update_yaxes(range=[0, 20],showgrid = True, gridwidth = 0.0002, gridcolor = 'grey')


#figN1.show()




# Plot!


st.subheader('All processing data view')

col1, col2 = st.columns(2)

with col1:
	st.plotly_chart(figrp, use_container_width=True)
with col2:
	st.plotly_chart(fig2b, use_container_width=True)


st.plotly_chart(figROX,use_container_width = True)
st.plotly_chart(figN1, use_container_width=True)
st.plotly_chart(fig1bbnbb, use_container_width=True)
# Streamlit widgets automatically run the script from top to bottom. Since
# this button is not connected to any other logic, it just causes a plain
# rerun.

def heat_map(df1, dye_choice, plate_choice):
    y = df1['Row_ID']
    x = df1['Col_ID']
    z = df1[dye_choice]
    fig = go.Figure(data=go.Heatmap(
        z=z,
        x=x,
        y=y,
        colorscale='magma'))
    fig.update_layout(title= (str(dye_choice) + ' HEATMAP: '+ str(plate_choice)),
                      xaxis_nticks=24,
                      yaxis_nticks = 16)
    fig.update_yaxes(autorange="reversed")
    st.plotly_chart(fig, use_container_width=True)
    

st.subheader('Individual array data - investigate data array by array')

plate = comp['Run_ID'].unique()
plate_choice = st.sidebar.selectbox('Select plate to analyse:', plate)
st.write(plate_choice)
select = comp[(comp.Run_ID == plate_choice)]
heat_map(select, 'ROX_RFU', plate_choice)
st.table(stats_ROX[stats_ROX['Run_ID'] == plate_choice])
#(stats_ROX[stats_ROX['Run_ID'] == plate_choice])
heat_map(select, 'norm_N_Cov', plate_choice)
st.table(stats_nFAM[stats_nFAM['Run_ID'] == plate_choice])
heat_map(select, 'norm_RNaseP', plate_choice)
st.table(stats_nCFO[stats_nCFO['Run_ID'] == plate_choice])



@st.cache
def convert_df(df):
     # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')

ROX = convert_df(stats_ROX)
FAM = convert_df(stats_FAM)
CFO = convert_df(stats_CFO)
nFAM = convert_df(stats_nFAM)
nCFO = convert_df(stats_nCFO)


st.sidebar.download_button(
     label="Download ROX data as CSV",
     data=ROX,
     file_name='ROX.csv',
     mime='text/csv',)

st.sidebar.download_button(
     label="Download FAM data as CSV",
     data=FAM,
     file_name='FAM.csv',
     mime='text/csv',)

st.sidebar.download_button(
     label="Download CFO data as CSV",
     data=CFO,
     file_name='CFO.csv',
     mime='text/csv',)

st.sidebar.download_button(
     label="Download nFAM data as CSV",
     data=nFAM,
     file_name='nFAM.csv',
     mime='text/csv',)

st.sidebar.download_button(
     label="Download nCFO data as CSV",
     data=nCFO,
     file_name='nCFO.csv',
     mime='text/csv',)
