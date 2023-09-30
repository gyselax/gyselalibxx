#!/usr/bin/env python3

# SPDX-License-Identifier: MIT

'''
file : plot_electric_field
date : 17/07/2023

Diagnostic for (X,Y,Vx,Vy) geometry

-----> Simple usage: python3 plot_electric_field.py ../../../build/simulations/geometryXYVxVy/landau/ --itime 10 
-----> Animation: python3 plot_electric_field.py ../../../build/simulations/geometryXYVxVy/landau/ --gif True

'''

# pylint: disable=invalid-name
# pylint: disable=wrong-import-order

from argparse import ArgumentParser
from pathlib import Path
from scipy import interpolate, stats

import matplotlib.pyplot as plt
import numpy as np
import os

from gysdata import DiskStore
from plot_utils import plot_field2d

import imageio.v2 as imageio

if __name__ == '__main__':
    parser = ArgumentParser(
        description='Plots the electric field')
    parser.add_argument('data_dir',
                        action='store',
                        nargs='?',
                        default=Path.cwd(),
                        type=Path,
                        help='location of the results')
    parser.add_argument('--gif',
                        action='store',
                        default=False,
                        type=bool,
                        help='produce animation')
    parser.add_argument('--itime',
                        action='store',
                        default=-1,
                        type=int,
                        help='time index')
    parser.add_argument('--grate',
                    metavar='<growth rate>',
                        action='store',
                        default=-0.1533,
                        type=float,
                        help='theoritical growth rate')

    args = parser.parse_args()

    # Load data
    theo_rate = args.grate
    data_structure_name = 'data_structure_XYVxVy.yaml'
    ds = DiskStore(args.data_dir, data_structure=Path(data_structure_name))
    epot = ds['electrostatic_potential']
    fdistribu = ds['fdistribu']
    Ex = -epot.differentiate('x')
    Ey = -epot.differentiate('y')
    Timer = epot.shape[0]
    time_init = epot.coords['time'].values[0]
    time_diag = epot.coords['time'].values[args.itime]
    time_final = epot.coords['time'].values[-1]

    output_filename = os.path.join(Path.cwd(), f'electric_energy_t{time_diag}.png')
    data_dict = {f't={time_init}': epot.sel(time=time_init),
                 f't={time_diag}': epot.sel(time=time_diag)}

    def plot_electric_potential(itime):
        '''
        Plot electric potential
        '''
        output_filename_epot = os.path.join(Path.cwd(), 'electric_potential.png')
        plot_field2d(epot[itime,:,:].T , '', output_filename_epot, cmap='jet', scale='linear')


    def plot_electric_field(itime):
        '''
        Plot electric field
        '''

        output_filename_Ex = os.path.join(Path.cwd(), 'electric_field_x.png')
        plot_field2d(Ex[itime-1,:,:].T , '', output_filename_Ex, cmap='jet', scale='linear')

        output_filename_Ey = os.path.join(Path.cwd(), 'electric_field_y.png')
        plot_field2d(Ey[itime-1,:,:].T , '', output_filename_Ey, cmap='jet', scale='linear')

    def plot_electric_energy_density_field(itime):
        '''
        Plot electric energy density field
        '''

        z = np.power(Ex[itime-1,:,:],2)+np.power(Ey[itime-1,:,:],2)
        output_filename_U = os.path.join(Path.cwd(), 'electric_energy_density_field.png')
        plot_field2d(z, '', output_filename_U, cmap='jet', scale='linear')


    def compute_electric_energy():
        '''
        Before displaying energy we have to perform a  two-dimensional integration.
        For that  we first interpolate the Electric field over the x,y plane.
        '''

        Energy = np.zeros(Timer)
        for k in range(1,Timer):
            z = np.power(Ex[k,:,:],2)+np.power(Ey[k,:,:],2)
            val = interpolate.RectBivariateSpline(epot.coords['x'].to_numpy(),
                epot.coords['y'].to_numpy(),z)
            Energy[k] = val.integral(epot.coords['x'].values[0], epot.coords['x'].values[-1],
                epot.coords['y'].values[0],epot.coords['y'].values[-1])
        return Energy

    def plot_electric_energy():
        '''
        Plot potential energy.
        '''

        Energy = compute_electric_energy()
        Emaxval = []
        tm = []
        for k in range(5,len(Energy)-2):
            if((Energy[k+1]-Energy[k-1])/Energy[k]< 0.05 and
                (Energy[k]-Energy[k-2])/Energy[k]> 0.1 and (Energy[k]-Energy[k+2])/Energy[k]>0.1):
                Emaxval.append(Energy[k])
                tm.append(epot.coords['time'].values[k])

        Emax = np.array(Emaxval)
        tm = np.array(tm)
        res = stats.linregress(tm,np.log(np.sqrt(Emax)))
        print(f'relative error on growth rate {(theo_rate-res[0]) /abs(res[0]):0.2f}')

        plt.clf()
        plt.plot(epot.coords['time'].values,np.log(np.sqrt(Energy)))
        plt.scatter(tm ,np.log(np.sqrt(Emax)),c='red')
        plt.plot(tm ,tm*res[0]+res[1]  ,c='green')
        plt.ylabel('Electric energy', fontsize=12)
        plt.xlabel('time', fontsize=12)
        plt.title(f'growthrate: {res[0]:0.3f} (theory: {theo_rate:0.3f})', fontsize=16)
        plt.savefig('Electric_energy.png')

    def plot_fdist(itime,Vx,Vy):
        '''
        Plot the initial distribution function in (X,Y) fo a fixed (Vx, Vy) velocity position.
        '''

        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot()# 111,projection='3d')
        im = ax.imshow(fdistribu[itime,0,:,:,Vx,Vy] )
        fig.colorbar(im)
        plt.title('initial distribution function at Vx='+str(Vx)+' Vy='+str(Vy))
        plt.savefig('test_fdistribution.png')

    if (not(args.gif)):
        plot_electric_potential(args.itime)
        plot_electric_field(args.itime)
        plot_electric_energy_density_field(args.itime)
        plot_electric_energy()
        plot_fdist(args.itime,31,31)

    else:
        min_efield = min([np.min(Ex),np.min(Ey)]).compute()
        max_efield = max([np.max(Ex),np.max(Ey)]).compute()
        min_epot = np.min(epot).compute()
        max_epot = np.max(epot).compute()
        min_energy_den = np.min(np.power(Ex,2)+np.power(Ey,2)).compute()
        max_energy_den = np.max(np.power(Ex,2)+np.power(Ey,2)).compute()

        images = []
        for itime in range(Timer):
            plot_field2d(epot[itime-1,:,:].T , '', os.path.join(Path.cwd(), 'temp.png'), vmin=min_epot, vmax=max_epot, cmap='jet', scale='linear')
            images.append(imageio.imread(os.path.join(Path.cwd(), 'temp.png')))
        imageio.mimsave('electric_potential.gif', images)

        images = []
        for itime in range(Timer):
            plot_field2d(Ex[itime-1,:,:].T , '', os.path.join(Path.cwd(), 'temp.png'), vmin=min_efield, vmax=max_efield, cmap='jet', scale='linear')
            images.append(imageio.imread(os.path.join(Path.cwd(), 'temp.png')))
        imageio.mimsave('electric_field_x.gif', images)

        images = []
        for itime in range(Timer):
            plot_field2d(Ey[itime-1,:,:].T , '', os.path.join(Path.cwd(), 'temp.png'), vmin=min_efield, vmax=max_efield, cmap='jet', scale='linear')
            images.append(imageio.imread(os.path.join(Path.cwd(), 'temp.png')))
        imageio.mimsave('electric_field_y.gif', images)

        images = []
        for itime in range(Timer):
            z = np.power(Ex[itime-1,:,:],2)+np.power(Ey[itime-1,:,:],2)
            plot_field2d(z, '', os.path.join(Path.cwd(), 'temp.png'), vmin=min_energy_den, vmax=max_energy_den, cmap='jet', scale='linear')
            images.append(imageio.imread(os.path.join(Path.cwd(), 'temp.png')))
        imageio.mimsave('electric_energy_density_field.gif', images)
