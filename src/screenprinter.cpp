// standard C++ library headers
#include <chrono>
#include <iostream>
#include <thread>

// the implemented class (last)
#include "screenprinter.hpp"

using std::cout;
using std::endl;
using std::flush;
using std::this_thread::sleep_for;
using std::chrono::milliseconds;

void ScreenPrinter::simulation_updated( const Distributed2DField& data )
{
	if ( data.distribution().rank() == 0 ) {
		cout << "at t="<<data.time()<<" : [" << endl;
	}
	cout<<flush; 
	MPI_Barrier(data.distribution().communicator());
	sleep_for(milliseconds(1));
	for ( int pyy = data.distribution().extent( DY )-1; pyy >=0 ; --pyy ) {
		for ( int yy = data.noghost_view().extent( DY )-1; yy >=0 ; --yy ) {
			for ( int pxx = 0; pxx < data.distribution().extent( DX ); ++pxx ) {
				if ( data.distribution().coord( DX ) == pxx && data.distribution().coord( DY ) == pyy ) {
					if ( pxx == 0 ) {
						cout << "  [";
					}
					for ( int xx = 0; xx < data.noghost_view().extent( DX ); ++xx ) {
						if ( 0 == data.noghost_view(yy, xx) ) {
							cout << " .";
						} else {
							cout << " " << data.noghost_view(yy, xx);
						}
					}
					if ( pxx == data.distribution().extent( DX ) - 1 ) {
						cout << " ]\n";
					}
				}
				cout<<""<<flush;
				MPI_Barrier(data.distribution().communicator());
				sleep_for(milliseconds(1));
			}
		}
	}
	if ( data.distribution().rank() == 0 ) {
		cout << "]" << endl;
	}
	cout<<flush;
	MPI_Barrier(data.distribution().communicator());
	sleep_for(milliseconds(1));
}
