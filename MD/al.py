import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib inline
a=4.046
center=[0,0,0]
global df
df = pd.DataFrame(columns=list("XYZ"))
# df.loc[len(df)]=[0,0,0]
def fcc(center):
  x1=center[0];y1=center[1];z1=center[2]
  x2=center[0]+a/2;y2=center[1]+a/2;z2=center[2]
  x3=center[0]+a/2;y3=center[1];z3=center[2]+a/2
  x4=center[0];y4=center[1]+a/2;z4=center[2]+a/2
  vx=0
  vy=0
  vz=0
  # atom=[x1,y1,z1]
  # df.loc[len(df)] = atom
  cell = pd.DataFrame([[x1,y1,z1,vx,vy,vz],[x2,y2,z2,vx,vy,vz],[x3,y3,z3,vx,vy,vz],[x4,y4,z4,vx,vy,vz]], columns=list('XYZxyz'))
 
  # df.append()
  # df.append()
  # df.append()
  # print(df3)
  return cell

def movecenter(nx,ny,nz):
  
  appended_data=[]
  for i in range(nx):
   for j in range(ny):
    for k in range(nz):
       center=[i*a,j*a,k*a]
       print(center)
       cell=fcc(center)
         #print(cell)
       appended_data.append(cell)
        #  pd.concat([df,df3],ignore_index=True)
  appended_data = pd.concat(appended_data)
  return appended_data

def dataview():
  df = pd.read_csv('data.csv', parse_dates=True)
  print(df.head())
  threedee = plt.figure().gca(projection='3d')
  threedee.scatter(df['X'], df['Y'], df['Z'])
  threedee.set_xlabel('X')
  threedee.set_ylabel('Y')
  threedee.set_zlabel('Z')
  plt.show()


appended_data=movecenter(3,3,3) 
appended_data=appended_data.drop_duplicates()
#appended_data=fcc(center)
# £print(t)
appended_data.to_csv('data.csv',index=False)
dataview()
