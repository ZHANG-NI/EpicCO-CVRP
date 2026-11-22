/*MIT License

Copyright(c) 2020 Thibaut Vidal

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files(the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions :

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/

#ifndef GENETIC_H
#define GENETIC_H

#include "Population.h"
#include "Individual.h"

#define CONSENSUS_K 5            // Elite size for the consensus don't-look marks (offer S5_v3)
#define CONSENSUS_REFRESH 100    // Iterations between refreshes of the consensus
#define DECOMP_MIN_CLIENTS 300   // Decomposition only for instances with at least this many customers
#define DECOMP_GROUP_SIZE 120    // Minimum number of customers in a subproblem
#define DECOMP_SUB_ITER 1000     // Non-improving iterations of the nested HGS
#define DECOMP_MIN_ITER 500      // Minimum number of main iterations between two decomposition phases
#define DECOMP_TIME_SHARE 0.3    // Initial share of the CPU time spent in decomposition
#define DECOMP_STALL 1000        // Decomposition only when the GA has not improved for this many iterations
#define DECOMP_SHARE_MIN 0.05
#define DECOMP_SHARE_MAX 0.6

class Genetic
{
public:

	Params & params;				// Problem parameters
	Split split;					// Split algorithm
	LocalSearch localSearch;		// Local Search structure
	Population population;			// Population (public for now to give access to the solutions, but should be be improved later on)
	Individual offspring;			// First individual to be used as input for the crossover

	double successEAX = 0.1;		// Smoothed rate at which EAX offspring improve on both parents (offer S6_v1)
	double successOX = 0.1;			// Smoothed rate at which OX offspring improve on both parents

	// EAX Crossover (returns false if no offspring could be produced) (offer S6_v1)
	bool crossoverEAX(Individual & result, const Individual & parentA, const Individual & parentB);

	// OX Crossover
	void crossoverOX(Individual & result, const Individual & parent1, const Individual & parent2);

	// SREX-like route exchange crossover: a block of polar-adjacent routes of parent2 replaces a block of routes of parent1
	void crossoverSREX(Individual & result, const Individual & parent1, const Individual & parent2);

	// Builds one SREX child (variant 0: duplicates removed from parent1's kept routes, variant 1: from parent2's moved routes)
	void buildSREXChild(Individual & child, const Individual & parent1, const Individual & parent2, const std::vector<int> & movedA, const std::vector<int> & movedB, int variant);

	// Working buffers for SREX
	// Inherited local optimality: tags the offspring routes identical to a route of parent A (1) or B (2) for the LS
	void tagInheritedRoutes(const Individual & child, const Individual & parentA, const Individual & parentB);
	std::vector<int> inhRoute, inhPos;
	std::vector<int> srexMarkA;		// 1 if customer belongs to the moved routes of parent1
	std::vector<int> srexMarkB;		// 1 if customer belongs to the moved routes of parent2
	std::vector<int> srexPlaced;	// 1 if customer already placed in the child
	Individual srexChild2;			// Second candidate child

	bool isSubproblem = false;		// Nested HGS run on a subproblem (offer S5_v1): stops after ap.nbIter non-improving iterations, no restart
	const Individual * initialSolution = nullptr; // Optional solution injected in the initial population

	// Route-based decomposition of the current best solution (groups of routes with close barycentre polar angles), each group solved by a nested HGS
	bool decompose();

    // Running the genetic algorithm until maxIterNonProd consecutive iterations or a time limit
    void run() ;

	// Constructor
	Genetic(Params & params);
};

#endif
