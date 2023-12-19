#!/usr/bin/env python
# coding: utf-8

# In[2]:


#loading Python Modules
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#Importing datatables
harris_p1 = pd.read_csv('C:/Users/Luke/Documents/University (UNSW)/Subjects/Physics/PHYS3116/PHYS3116 Assessment Tasks/Computational Assignment/HarrisPartI.csv')
harris_p3 = pd.read_csv('C:/Users/Luke/Documents/University (UNSW)/Subjects/Physics/PHYS3116/PHYS3116 Assessment Tasks/Computational Assignment/HarrisPartIII.csv')
gaia_df = pd.read_csv('C:/Users/Luke/Documents/University (UNSW)/Subjects/Physics/PHYS3116/PHYS3116 Assessment Tasks/Computational Assignment/result.csv')


# In[3]:


#Function to swap values in two columns based on alphabetical order
def swap_columns_alphabetically(df, col1, col2):
    for index, row in df.iterrows():
        if row[col1]==type(str) or row[col2]==type(str):
             row[col1].tostr
        if row[col1] > row[col2]:
            df.at[index, col1], df.at[index, col2] = row[col2], row[col1]


# In[4]:


#Function to find distance in X for stars from galactic centre
def custom_transform(x):
    if x < 8:
        return 8 - x
    else:
        return x + 8


# In[5]:


def angle_transform(x):
    if x> np.pi/2:
        return x-(np.pi/2)
    else:
        return (np.pi/2)-x


# In[6]:


def angle_transform2(x):
    if x<np.pi/2:
        return x
    else:
        return x-(np.pi/2)


# In[7]:


#Merge Harris_PartI and Harris_PartIII
harris_df = pd.merge(harris_p1, harris_p3, on='ID')

#Replace nan in harris data to empty string
harris_df['Name']=harris_df['Name'].fillna('')
#harris_df.to_csv('harris_dffillna.csv', index=False)


# In[8]:


# Drop NaN values for v_r as we need recessional velocity to have a value for calculations
harris_df= harris_df.dropna(subset='v_r')


# In[9]:


#Changing header names in gaia_DF to 'ID' and 'Name'
new_header = ['ID','Name','RAJ2000','DEJ2000','pmRA','e_pmRA','pmDE','e_pmDE','corr','plx','e_plx','Rscale','Nstar','SimbadName']
gaia_df.columns=new_header

# Call the function to alphabetise columns 'ID' and 'Name' of the DataFrames
swap_columns_alphabetically(gaia_df, 'ID', 'Name')
swap_columns_alphabetically(harris_df, 'ID', 'Name')


# In[10]:


# Merge harris_df and gaia_df
MW_df = pd.merge(gaia_df, harris_df, on='Name')

# Delete duplicate columns with secondary name 'ID_x' and 'ID_y'
columns_to_delete = ['ID_x', 'ID_y']
MW_df = MW_df.drop(columns=columns_to_delete)


# In[11]:


#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
#xx MODEL 1: Total Velocity as a plot of galactocentric distance xx
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# total value of proper motion
MW_df['pm_total'] = np.sqrt( (MW_df['pmRA'])**2+(MW_df['pmDE'])**2)

# linear velocity for proper motion in km/s
MW_df['v_pm']=(4.74*MW_df['pm_total']*(MW_df['R_Sun']))

# Find the total magnitude of velocity by adding the recessional vel and tangential vel
MW_df['v']= np.sqrt( (MW_df['v_pm']**2) + (MW_df['v_r']**2))


# In[13]:


plt.scatter(MW_df['R_gc'], MW_df['v'])
plt.axvline(x=40, color='r', ls='--')
plt.ylabel('Total Velocity (km/s)')
plt.xlabel('Radius from galactic centre (kpc)')
plt.title('Scatter Plot of Total Vel. and Radius to Galactic Centre')
plt.show()


# In[14]:


plt.scatter(MW_df['R_gc'], MW_df['v'])
plt.xlim(0, 40)
plt.ylabel('Total Velocity (km/s)')
plt.xlabel('Radius from galactic centre (kpc)')
plt.title('Scatter Plot of Total Vel. and Radius to Galactic Centre')
plt.show()


# In[15]:


mean_total_vel=np.mean(MW_df['v'])
velocity_dispersion = np.std(MW_df['v'], ddof=1)

print(f"Mean Velocity: {mean_total_vel}")
print(f"Velocity Dispersion: {velocity_dispersion}")


# In[16]:


#Find distance from sun to globular clusters and distance from galactic centre(gc) to globular clusters
MW_df['R_sun_xy'] = np.sqrt(MW_df['X']**2 + MW_df['Y']**2)

# Find stars altitude relative to sun and vector component of recessional velocity (v_r) in XY-Plane
MW_df['alt_sun'] = np.arctan( MW_df['Z'] / MW_df['R_sun_xy'])
MW_df['v_r_xy'] = MW_df['v_r']*np.cos(MW_df['alt_sun'])

MW_df['X_gc'] = MW_df['X'].apply(custom_transform) # is x < than or >= to 8kpc: inside or outside Solar orbit.
MW_df['DCB']=np.arctan( MW_df['Y']/MW_df['X'])
MW_df['DAB']=np.arctan( MW_df['Y']/MW_df['X_gc'])
MW_df['ABC']=(abs(MW_df['DCB']) + abs(MW_df['DAB']))
MW_df['phi']=MW_df['ABC'].apply(angle_transform)
#Recessional Velocity component tangential to galaxy centre
MW_df['V_v_r']=np.sin(MW_df['phi'])*MW_df['v_r_xy']


# In[17]:


# Velocity of proper motion right ascension in km/s
MW_df['v_pmra']=4.74*MW_df['R_gc']*MW_df['pmRA']

MW_df['theta']=MW_df['ABC'].apply(angle_transform2)
#Proper motion velocity component tangential to galaxy centre
MW_df['V_pm']=np.sin(MW_df['theta'])*MW_df['v_pmra']

#Total Tangential Velocity
MW_df['velocity_V']=MW_df['V_pm']+MW_df['V_v_r']


# In[18]:


plt.scatter(MW_df['R_gc'], MW_df['velocity_V'])
plt.xlim(0, 40)
plt.ylabel('Total Velocity (km/s)')
plt.xlabel('Radius from galactic centre (kpc)')
plt.title('Scatter Plot of Total Vel. and Radius to Galactic Centre')
plt.show()


# In[19]:


#Bulk Rotation Check
same_sign = 0
different_sign = 0
dummy_variable = 0

for index, row in MW_df.head(200).iterrows():
    if row['Y'] >0 and row['v_r'] > 0:
        same_sign+=1
    else:
        different_sign +=1

print("\"Clockwise\" count", same_sign)
print("\"Anti-clockwise\" count:", different_sign)
print("Percentage:", same_sign/different_sign*100,"%.")


# In[ ]:




