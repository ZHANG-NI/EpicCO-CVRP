#include "Genetic.h"

void Genetic::run()
{	
	/* INITIAL POPULATION */
	population.generatePopulation();
	if (initialSolution != nullptr) population.addIndividual(*initialSolution, false);
	// Adaptive control of the decomposition (offer S5_v1): its share of the CPU time follows its recent rate of improvement relative to the main search
	double decompShare = DECOMP_TIME_SHARE;
	double lastPhaseDuration = 0.;   // CPU time of the last decomposition phase
	double lastPhaseEnd = 0.;        // Time stamp of the end of the last decomposition phase
	double costAtLastPhaseEnd = 1.e30;
	int lastDecompIter = 0;

	int nbIter;
	int nbIterNonProd = 1;
	double phaseStart = 0.;     // CPU time at which the current HGS phase (population) started
	double phaseBestTime = 0.;  // CPU time of the last improvement of the current phase
	if (params.verbose) std::cout << "----- STARTING GENETIC ALGORITHM" << std::endl;
	for (nbIter = 0 ; (nbIterNonProd <= params.ap.nbIter || (!isSubproblem && params.ap.timeLimit != 0)) && (params.ap.timeLimit == 0 || (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC < params.ap.timeLimit) ; nbIter++)
	{	
		/* ELITE CONSENSUS: refresh the don't-look marks of the local search (offer S5_v3) */
		if (nbIter % CONSENSUS_REFRESH == 0) population.computeConsensus(CONSENSUS_K, localSearch.consensusA, localSearch.consensusB);

		/* SELECTION AND CROSSOVER */
		const Individual & parentA = population.getBinaryTournament();
		const Individual & parentB = population.getBinaryTournament();
		const double parentsCost = std::min<double>(parentA.eval.penalizedCost, parentB.eval.penalizedCost);
		int usedOp = 0; // 0: OX, 1: EAX, 2: SREX
		if (params.ran() % 100 < (unsigned)params.ap.srexRate) { crossoverSREX(offspring, parentA, parentB); usedOp = 2; }
		else
		{
			double probEAX = std::min<double>(0.9, std::max<double>(0.1, successEAX / (successEAX + successOX)));
			bool usedEAX = (std::uniform_real_distribution<double>(0., 1.)(params.ran) < probEAX && &parentA != &parentB && crossoverEAX(offspring, parentA, parentB));
			if (!usedEAX) crossoverOX(offspring, parentA, parentB);
			usedOp = usedEAX ? 1 : 0;
		}

		/* LOCAL SEARCH */
		tagInheritedRoutes(offspring, parentA, parentB);
		localSearch.run(offspring, params.penaltyCapacity, params.penaltyDuration);
		/* ADAPTIVE EAX/OX SELECTION (offer S6_v1): exponentially smoothed rate of offspring better than both parents after education */
		{
			double reward = (offspring.eval.penalizedCost < parentsCost - MY_EPSILON) ? 1.0 : 0.0;
			if (usedOp == 1) successEAX = 0.995 * successEAX + 0.005 * reward;
			else if (usedOp == 0) successOX = 0.995 * successOX + 0.005 * reward;
		}
		bool isNewBest = population.addIndividual(offspring,true);
		if (!offspring.eval.isFeasible && params.ran()%2 == 0) // Repair half of the solutions in case of infeasibility
		{
			// The offspring is a local optimum under the lower penalty: moves among its feasible routes cannot become
			// improving under a higher penalty, so the first LS loop only evaluates moves involving an infeasible route
			for (int r = 0; r < params.nbVehicles; r++)
			{
				double load = 0.;
				for (int c : offspring.chromR[r]) load += params.cli[c].demand;
				localSearch.routeTag[r] = (load <= params.vehicleCapacity + MY_EPSILON && !params.isDurationConstraint) ? 1 : 0;
			}
			localSearch.hasRouteTags = true;
			localSearch.run(offspring, params.penaltyCapacity*10., params.penaltyDuration*10.);
			if (offspring.eval.isFeasible) isNewBest = (population.addIndividual(offspring,false) || isNewBest);
		}

		/* TRACKING THE NUMBER OF ITERATIONS SINCE LAST SOLUTION IMPROVEMENT */
		if (isNewBest) { nbIterNonProd = 1; phaseBestTime = (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC; }
		else nbIterNonProd ++ ;

		/* DIVERSIFICATION, PENALTY MANAGEMENT AND TRACES */
		if (nbIter % params.ap.nbIterPenaltyManagement == 0) population.managePenalties();
		if (nbIter % params.ap.nbIterTraces == 0) population.printState(nbIter, nbIterNonProd);

		/* FOR TESTS INVOLVING SUCCESSIVE RUNS UNTIL A TIME LIMIT: WE RESET THE ALGORITHM/POPULATION EACH TIME maxIterNonProd IS ATTAINED*/
		/* ROUTE-BASED DECOMPOSITION, with a CPU time share adapted to its success */
		if (!isSubproblem && params.ap.timeLimit != 0 && params.nbClients >= DECOMP_MIN_CLIENTS && nbIter - lastDecompIter >= DECOMP_MIN_ITER && nbIterNonProd >= DECOMP_STALL
			&& population.getBestFeasible() != NULL)
		{
			double now = (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC;
			if (now - lastPhaseEnd >= lastPhaseDuration * (1. - decompShare) / decompShare)
			{
				double costBefore = population.getBestFeasible()->eval.penalizedCost;
				double mainRate = (costAtLastPhaseEnd < 1.e29) ? std::max<double>(0., costAtLastPhaseEnd - costBefore) / std::max<double>(1.e-3, now - lastPhaseEnd) : -1.;
				if (decompose()) { nbIterNonProd = 1; phaseBestTime = (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC; }
				double end = (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC;
				lastPhaseDuration = end - now;
				lastPhaseEnd = end;
				costAtLastPhaseEnd = population.getBestFeasible()->eval.penalizedCost;
				double decompRate = std::max<double>(0., costBefore - costAtLastPhaseEnd) / std::max<double>(1.e-3, lastPhaseDuration);
				if (mainRate >= 0.)
				{
					double ratio = (decompRate + mainRate > 0.) ? decompRate / (decompRate + mainRate) : decompShare * 0.5;
					decompShare = std::min<double>(DECOMP_SHARE_MAX, std::max<double>(DECOMP_SHARE_MIN, 0.5 * decompShare + 0.5 * ratio));
				}
				lastDecompIter = nbIter;
				if (params.verbose) std::printf("DECOMP-CTRL | rates main %.3f decomp %.3f | share %.3f\n", mainRate, decompRate, decompShare);
			}
		}

		// Restart after nbIter/2 non-improving iterations (offer S6_v2), but only if the remaining time is at least the time
		// the current phase needed to reach its best (otherwise a new population cannot catch up: keep the current one)
		if (!isSubproblem && params.ap.timeLimit != 0 && nbIterNonProd == params.ap.nbIter / 2
			&& params.ap.timeLimit - (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC >= phaseBestTime - phaseStart)
		{
			phaseStart = phaseBestTime = (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC;
			std::fill(localSearch.consensusA.begin(), localSearch.consensusA.end(), -1); // no stale marks while the new population is built
			population.restart();
			population.computeConsensus(CONSENSUS_K, localSearch.consensusA, localSearch.consensusB);
			nbIterNonProd = 1;
		}
	}
	if (params.verbose) std::cout << "----- GENETIC ALGORITHM FINISHED AFTER " << nbIter << " ITERATIONS. TIME SPENT: " << (double)(clock() - params.startTime) / (double)CLOCKS_PER_SEC << std::endl;
}

void Genetic::tagInheritedRoutes(const Individual & child, const Individual & parentA, const Individual & parentB)
{
	std::vector<int> & tag = localSearch.routeTag;
	std::fill(tag.begin(), tag.end(), 0);
	if ((int)inhRoute.size() != params.nbClients + 1) { inhRoute.assign(params.nbClients + 1, 0); inhPos.assign(params.nbClients + 1, 0); }
	for (int k = 1; k <= 2; k++)
	{
		const Individual & P = (k == 1) ? parentA : parentB;
		if (k == 2 && &parentB == &parentA) break;
		for (int r = 0; r < params.nbVehicles; r++)
			for (int i = 0; i < (int)P.chromR[r].size(); i++) { inhRoute[P.chromR[r][i]] = r; inhPos[P.chromR[r][i]] = i; }
		for (int r = 0; r < params.nbVehicles; r++)
		{
			const std::vector<int> & R = child.chromR[r];
			if (tag[r] != 0 || R.empty()) continue;
			const std::vector<int> & PR = P.chromR[inhRoute[R[0]]];
			int s = (int)R.size();
			if ((int)PR.size() != s) continue;
			bool same = false;
			if (inhPos[R[0]] == 0)
			{
				same = true;
				for (int i = 1; i < s && same; i++) same = (PR[i] == R[i]);
			}
			if (!same && params.isSymmetric && inhPos[R[0]] == s - 1)
			{
				same = true;
				for (int i = 1; i < s && same; i++) same = (PR[s - 1 - i] == R[i]);
			}
			if (same) tag[r] = k;
		}
	}
	localSearch.hasRouteTags = true;
}

bool Genetic::decompose()
{
	const Individual * best = population.getBestFeasible();
	if (best == NULL) return false;

	// Barycentres of the non-empty routes
	std::vector<int> rid;
	std::vector<double> bxs, bys;
	for (int r = 0; r < params.nbVehicles; r++)
	{
		if (best->chromR[r].empty()) continue;
		double bx = 0., by = 0.;
		for (int c : best->chromR[r]) { bx += params.cli[c].coordX; by += params.cli[c].coordY; }
		rid.push_back(r);
		bxs.push_back(bx / (double)best->chromR[r].size());
		bys.push_back(by / (double)best->chromR[r].size());
	}
	int nbR = (int)rid.size();
	if (nbR < 2) return false;

	std::vector<std::vector<int>> groups;
	if (params.ran() % 2 == 0)
	{
		// Polar wedges: routes ordered by the polar angle of their barycentre, consecutive routes grouped until each
		// group has at least DECOMP_GROUP_SIZE customers, from a random starting route
		std::vector<std::pair<double,int>> angles;
		for (int k = 0; k < nbR; k++) angles.push_back({atan2(bys[k] - params.cli[0].coordY, bxs[k] - params.cli[0].coordX), rid[k]});
		std::sort(angles.begin(), angles.end());
		int offset = params.ran() % nbR;
		std::vector<int> current;
		int count = 0;
		for (int k = 0; k < nbR; k++)
		{
			int r = angles[(k + offset) % nbR].second;
			current.push_back(r);
			count += (int)best->chromR[r].size();
			if (count >= DECOMP_GROUP_SIZE) { groups.push_back(current); current.clear(); count = 0; }
		}
		if (!current.empty())
		{
			if (!groups.empty() && count < DECOMP_GROUP_SIZE / 2) groups.back().insert(groups.back().end(), current.begin(), current.end());
			else groups.push_back(current);
		}
	}
	else
	{
		// Compact clusters: from a random unassigned seed route, add the unassigned routes with the closest barycentres
		// until the group has at least DECOMP_GROUP_SIZE customers (remaining routes are absorbed if they are too few)
		std::vector<bool> assigned(nbR, false);
		int remainingCustomers = 0;
		for (int k = 0; k < nbR; k++) remainingCustomers += (int)best->chromR[rid[k]].size();
		int nbAssigned = 0;
		while (nbAssigned < nbR)
		{
			std::vector<int> unassigned;
			for (int k = 0; k < nbR; k++) if (!assigned[k]) unassigned.push_back(k);
			int seed = unassigned[params.ran() % unassigned.size()];
			std::vector<std::pair<double,int>> byDist;
			for (int k : unassigned) byDist.push_back({(bxs[k]-bxs[seed])*(bxs[k]-bxs[seed]) + (bys[k]-bys[seed])*(bys[k]-bys[seed]), k});
			std::sort(byDist.begin(), byDist.end());
			std::vector<int> current;
			int count = 0;
			for (auto & e : byDist)
			{
				if (count >= DECOMP_GROUP_SIZE && remainingCustomers - count >= DECOMP_GROUP_SIZE / 2) break;
				current.push_back(rid[e.second]);
				assigned[e.second] = true;
				nbAssigned++;
				count += (int)best->chromR[rid[e.second]].size();
			}
			remainingCustomers -= count;
			groups.push_back(current);
		}
	}
	if (groups.size() < 2) return false;

	Individual newSol = *best;
	newSol.indivsPerProximity.clear();
	bool improved = false;
	for (const std::vector<int> & group : groups)
	{
		double remaining = params.ap.timeLimit - (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC;
		if (remaining < 0.5) break;

		// Building the subproblem (local index 0 = depot)
		std::vector<int> globalId = {0};
		double subDemand = 0.;
		for (int r : group) for (int c : newSol.chromR[r]) { globalId.push_back(c); subDemand += params.cli[c].demand; }
		int m = (int)globalId.size() - 1;
		std::vector<double> x(m+1), y(m+1), serv(m+1), dem(m+1);
		std::vector<std::vector<double>> dist(m+1, std::vector<double>(m+1));
		for (int i = 0; i <= m; i++)
		{
			const Client & cl = params.cli[globalId[i]];
			x[i] = cl.coordX; y[i] = cl.coordY; serv[i] = cl.serviceDuration; dem[i] = cl.demand;
			for (int j = 0; j <= m; j++) dist[i][j] = params.timeCost[globalId[i]][globalId[j]];
		}
		AlgorithmParameters sap = params.ap;
		sap.nbIter = DECOMP_SUB_ITER;
		sap.seed = (int)(params.ran() % 1000000007);
		sap.timeLimit = std::max<double>(0.1, remaining - 0.3);
		int subVeh = std::max<int>((int)std::ceil(1.3*subDemand/params.vehicleCapacity) + 3, (int)group.size());
		Params sp(x, y, dist, serv, dem, params.vehicleCapacity, params.durationLimit, subVeh, params.isDurationConstraint, false, sap);

		// Current routes of the group, as the initial solution of the subproblem
		Individual init(sp);
		int pos = 0, k = 0, local = 1;
		for (int r : group)
		{
			for (int idx = 0; idx < (int)newSol.chromR[r].size(); idx++) { init.chromR[k].push_back(local); init.chromT[pos++] = local; local++; }
			k++;
		}
		init.evaluateCompleteCost(sp);

		Genetic sub(sp);
		sub.isSubproblem = true;
		sub.initialSolution = &init;
		sub.run();
		const Individual * sb = sub.population.getBestFound();
		if (sb == NULL || sb->eval.penalizedCost > init.eval.penalizedCost - MY_EPSILON) continue;

		// Writing back the improved routes in the route slots of the group (and empty slots if needed)
		std::vector<int> slots = group;
		for (int r : group) newSol.chromR[r].clear();
		for (int r = 0; r < params.nbVehicles; r++) if (newSol.chromR[r].empty() && std::find(group.begin(), group.end(), r) == group.end()) slots.push_back(r);
		int s = 0;
		for (int r = 0; r < sp.nbVehicles; r++)
		{
			if (sb->chromR[r].empty()) continue;
			if (s >= (int)slots.size()) throw std::string("Decomposition: not enough route slots");
			for (int c : sb->chromR[r]) newSol.chromR[slots[s]].push_back(globalId[c]);
			s++;
		}
		improved = true;
	}
	if (!improved) return false;

	int pos = 0;
	for (int r = 0; r < params.nbVehicles; r++) for (int c : newSol.chromR[r]) newSol.chromT[pos++] = c;
	newSol.evaluateCompleteCost(params);
	double before = newSol.eval.penalizedCost;
	double bestCost = best->eval.penalizedCost;
	bool res = population.addIndividual(newSol, false); // The recombined solution is feasible: keep it as is
	localSearch.run(newSol, params.penaltyCapacity, params.penaltyDuration); // and also its improvement across group boundaries
	if (newSol.eval.penalizedCost < before - MY_EPSILON) res = population.addIndividual(newSol, false) || res;
	if (params.verbose) std::printf("DECOMP | groups %d | best %.2f | recombined %.2f | after LS %.2f (feas %d) | T(s) %.2f\n", (int)groups.size(), bestCost, before, newSol.eval.penalizedCost, (int)newSol.eval.isFeasible, (double)(clock()-params.startTime)/(double)CLOCKS_PER_SEC);
	return res;
}

void Genetic::crossoverOX(Individual & result, const Individual & parent1, const Individual & parent2)
{
	// Frequency table to track the customers which have been already inserted
	std::vector <bool> freqClient = std::vector <bool> (params.nbClients + 1, false);

	// Picking the beginning and end of the crossover zone
	std::uniform_int_distribution<> distr(0, params.nbClients-1);
	int start = distr(params.ran);
	int end = distr(params.ran);

	// Avoid that start and end coincide by accident
	while (end == start) end = distr(params.ran);

	// Copy from start to end
	int j = start;
	while (j % params.nbClients != (end + 1) % params.nbClients)
	{
		result.chromT[j % params.nbClients] = parent1.chromT[j % params.nbClients];
		freqClient[result.chromT[j % params.nbClients]] = true;
		j++;
	}

	// Fill the remaining elements in the order given by the second parent
	for (int i = 1; i <= params.nbClients; i++)
	{
		int temp = parent2.chromT[(end + i) % params.nbClients];
		if (freqClient[temp] == false)
		{
			result.chromT[j % params.nbClients] = temp;
			j++;
		}
	}

	// Complete the individual with the Split algorithm
	split.generalSplit(result, parent1.eval.nbRoutes);
}

static int nbNonEmptyRoutes(const Individual & indiv)
{
	int nb = 0;
	for (const auto & r : indiv.chromR) if (!r.empty()) nb++;
	return nb;
}

void Genetic::buildSREXChild(Individual & child, const Individual & parent1, const Individual & parent2, const std::vector<int> & movedA, const std::vector<int> & movedB, int variant)
{
	int nA = nbNonEmptyRoutes(parent1);
	std::fill(srexPlaced.begin(), srexPlaced.end(), 0);
	for (auto & r : child.chromR) r.clear();
	int k = 0;

	// Routes of parent1 which are kept
	for (int i = 0; i < nA; i++)
	{
		bool moved = false;
		for (int r : movedA) if (r == i) moved = true;
		if (moved) continue;
		for (int c : parent1.chromR[i])
			if (variant == 1 || !srexMarkB[c]) { child.chromR[k].push_back(c); srexPlaced[c] = 1; }
		if (!child.chromR[k].empty()) k++;
	}
	// Routes of parent2 which are moved in
	for (int rB : movedB)
	{
		for (int c : parent2.chromR[rB])
			if (!srexPlaced[c]) { child.chromR[k].push_back(c); srexPlaced[c] = 1; }
		if (!child.chromR[k].empty()) k++;
	}

	// Load of each route
	std::vector<double> load(k, 0.);
	for (int r = 0; r < k; r++) for (int c : child.chromR[r]) load[r] += params.cli[c].demand;

	// Greedy best insertion (distance + capacity penalty) of the customers which are still missing
	for (int c = 1; c <= params.nbClients; c++)
	{
		if (srexPlaced[c]) continue;
		double bestCost = 1.e30; int bestR = -1; int bestPos = -1;
		for (int r = 0; r < k; r++)
		{
			double penaltyDelta = params.penaltyCapacity * (std::max<double>(0., load[r] + params.cli[c].demand - params.vehicleCapacity) - std::max<double>(0., load[r] - params.vehicleCapacity));
			if (penaltyDelta >= bestCost) continue;
			const std::vector<int> & route = child.chromR[r];
			int prev = 0;
			for (int pos = 0; pos <= (int)route.size(); pos++)
			{
				int next = (pos < (int)route.size()) ? route[pos] : 0;
				double cost = params.timeCost[prev][c] + params.timeCost[c][next] - params.timeCost[prev][next] + penaltyDelta;
				if (cost < bestCost) { bestCost = cost; bestR = r; bestPos = pos; }
				prev = next;
			}
		}
		// Possibility to open a new route if some vehicles are available
		if (k < params.nbVehicles)
		{
			double cost = params.timeCost[0][c] + params.timeCost[c][0];
			if (cost < bestCost) { bestCost = cost; bestR = k; bestPos = 0; k++; load.push_back(0.); }
		}
		child.chromR[bestR].insert(child.chromR[bestR].begin() + bestPos, c);
		load[bestR] += params.cli[c].demand;
		srexPlaced[c] = 1;
	}

	// Giant tour consistent with the routes
	int pos = 0;
	for (int r = 0; r < params.nbVehicles; r++)
		for (int c : child.chromR[r]) child.chromT[pos++] = c;
	child.evaluateCompleteCost(params);
}

void Genetic::crossoverSREX(Individual & result, const Individual & parent1, const Individual & parent2)
{
	int nA = nbNonEmptyRoutes(parent1);
	int nB = nbNonEmptyRoutes(parent2);
	if (nA < 2 || nB < 2) { crossoverOX(result, parent1, parent2); return; }

	// Routes are stored by increasing polar angle of their barycenter (see LocalSearch::exportIndividual), so consecutive indices are spatial neighbors
	int startA = params.ran() % nA;
	int nbMoved = 1 + params.ran() % std::max<int>(1, std::min<int>(nA, nB) / 2);
	int startB = params.ran() % nB;
	std::vector<int> movedA, movedB;
	if (params.ran() % 2 == 0)
	{
		// Compact blocks: the seed route of parent1 and its nearest routes by barycentre; the routes of parent2 with the
		// best overlap (customers inside minus customers outside parent1's block)
		std::vector<double> bx(nA, 0.), by(nA, 0.);
		for (int r = 0; r < nA; r++)
		{
			for (int c : parent1.chromR[r]) { bx[r] += params.cli[c].coordX; by[r] += params.cli[c].coordY; }
			bx[r] /= parent1.chromR[r].size(); by[r] /= parent1.chromR[r].size();
		}
		std::vector<std::pair<double,int>> byDist;
		for (int r = 0; r < nA; r++) byDist.push_back({(bx[r]-bx[startA])*(bx[r]-bx[startA]) + (by[r]-by[startA])*(by[r]-by[startA]) + (r == startA ? -1. : 0.), r});
		std::partial_sort(byDist.begin(), byDist.begin() + nbMoved, byDist.end());
		std::fill(srexMarkA.begin(), srexMarkA.end(), 0);
		for (int j = 0; j < nbMoved; j++) { movedA.push_back(byDist[j].second); for (int c : parent1.chromR[byDist[j].second]) srexMarkA[c] = 1; }
		std::vector<std::pair<int,int>> scoreB;
		for (int r = 0; r < nB; r++)
		{
			int sc = 0;
			for (int c : parent2.chromR[r]) sc += srexMarkA[c] ? 1 : -1;
			scoreB.push_back({-sc, r});
		}
		std::partial_sort(scoreB.begin(), scoreB.begin() + nbMoved, scoreB.end());
		std::fill(srexMarkB.begin(), srexMarkB.end(), 0);
		for (int j = 0; j < nbMoved; j++) { movedB.push_back(scoreB[j].second); for (int c : parent2.chromR[scoreB[j].second]) srexMarkB[c] = 1; }
		buildSREXChild(result, parent1, parent2, movedA, movedB, 0);
		buildSREXChild(srexChild2, parent1, parent2, movedA, movedB, 1);
		if (srexChild2.eval.penalizedCost < result.eval.penalizedCost - MY_EPSILON) result = srexChild2;
		return;
	}

	// Choose startB such that the moved block of parent2 covers the customers of the moved block of parent1 as well as possible
	std::fill(srexMarkA.begin(), srexMarkA.end(), 0);
	for (int j = 0; j < nbMoved; j++) for (int c : parent1.chromR[(startA + j) % nA]) srexMarkA[c] = 1;
	auto overlapScore = [&](int sB) {
		int score = 0;
		for (int j = 0; j < nbMoved; j++) for (int c : parent2.chromR[(sB + j) % nB]) score += srexMarkA[c] ? 1 : -1;
		return score;
	};
	// Initial guess: the route of parent2 containing the most customers of the first moved route of parent1
	{
		std::vector<int> cnt(nB, 0);
		std::vector<int> routeOf(params.nbClients + 1, -1);
		for (int r = 0; r < nB; r++) for (int c : parent2.chromR[r]) routeOf[c] = r;
		for (int c : parent1.chromR[startA]) cnt[routeOf[c]]++;
		startB = (int)(std::max_element(cnt.begin(), cnt.end()) - cnt.begin());
	}
	int bestScore = overlapScore(startB);
	bool improved = true;
	while (improved)
	{
		improved = false;
		for (int d : {-1, 1})
		{
			int sB = (startB + d + nB) % nB;
			int sc = overlapScore(sB);
			if (sc > bestScore) { bestScore = sc; startB = sB; improved = true; break; }
		}
	}
	std::fill(srexMarkB.begin(), srexMarkB.end(), 0);
	for (int j = 0; j < nbMoved; j++) for (int c : parent2.chromR[(startB + j) % nB]) srexMarkB[c] = 1;

	for (int j = 0; j < nbMoved; j++) { movedA.push_back((startA + j) % nA); movedB.push_back((startB + j) % nB); }
	buildSREXChild(result, parent1, parent2, movedA, movedB, 0);
	buildSREXChild(srexChild2, parent1, parent2, movedA, movedB, 1);
	if (srexChild2.eval.penalizedCost < result.eval.penalizedCost - MY_EPSILON) result = srexChild2;
}

Genetic::Genetic(Params & params) : 
	params(params), 
	split(params),
	localSearch(params),
	population(params,this->split,this->localSearch),
	offspring(params),
	srexMarkA(params.nbClients + 1),
	srexMarkB(params.nbClients + 1),
	srexPlaced(params.nbClients + 1),
	srexChild2(params){}

