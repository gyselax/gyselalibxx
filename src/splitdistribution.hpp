#pragma once

#include <mpi.h>

#include <array>

class BlockDistribution
{
public:
	/// number of blocks in this dimension
	int size() const;
	
	/// The communicator parallel to this distribution
	MPI_Comm m_parallel_comm;
	
	/// The communicator orthogonal to this distribution
	MPI_Comm m_orth_comm;
	
};

template<int N>
class SplitDistribution
{
	std::array<BlockDistribution, N> m_distributions;
	
	/// The comm with all processes taking part in the distribution
	MPI_Comm m_comm;
	
public:
	
	
};
