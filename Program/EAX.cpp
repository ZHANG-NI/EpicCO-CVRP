// Edge assembly crossover (EAX) for the CVRP, used as an alternative to the OX crossover.
// AB-cycles are built from the symmetric difference of the edge multisets of two parents (the depot has degree 2m);
// a spatially localized block of AB-cycles is applied to parent A; the resulting subtours (cycles without the depot)
// are merged into routes by a 2-opt-style reconnection, and the offspring is then educated by the local search.
#include "Genetic.h"

namespace {
struct ABEdge { int u, v; bool isA; };

void removeOne(std::vector<int> & lst, int x)
{
	for (int i = 0; i < (int)lst.size(); i++)
		if (lst[i] == x) { lst[i] = lst.back(); lst.pop_back(); return; }
}
}

bool Genetic::crossoverEAX(Individual & result, const Individual & parentA, const Individual & parentB)
{
	const int N = params.nbClients + 1;
	std::vector<std::vector<int>> adjA(N), adjB(N);
	auto buildAdj = [&](const Individual & I, std::vector<std::vector<int>> & adj)
	{
		for (const auto & route : I.chromR)
		{
			if (route.empty()) continue;
			int prev = 0;
			for (int c : route) { adj[prev].push_back(c); adj[c].push_back(prev); prev = c; }
			adj[prev].push_back(0); adj[0].push_back(prev);
		}
	};
	buildAdj(parentA, adjA);
	buildAdj(parentB, adjB);

	// Symmetric difference (with multiplicities)
	std::vector<std::vector<int>> dA(N), dB(N);
	std::vector<int> cnt(N, 0);
	int nbDiff = 0;
	for (int v = 0; v < N; v++)
	{
		for (int w : adjB[v]) cnt[w]++;
		for (int w : adjA[v]) { if (cnt[w] > 0) cnt[w]--; else dA[v].push_back(w); }
		for (int w : adjB[v]) cnt[w] = 0;
		for (int w : adjA[v]) cnt[w]++;
		for (int w : adjB[v]) { if (cnt[w] > 0) cnt[w]--; else dB[v].push_back(w); }
		for (int w : adjA[v]) cnt[w] = 0;
		nbDiff += (int)dA[v].size();
	}
	if (nbDiff == 0) return false;

	// Decomposition of the difference into AB-cycles by random alternating walks
	std::vector<std::vector<ABEdge>> cycles;
	std::vector<int> order(N - 1);
	for (int i = 0; i < N - 1; i++) order[i] = i + 1;
	std::shuffle(order.begin(), order.end(), params.ran);
	std::vector<int> P;
	std::vector<char> T;
	std::vector<std::vector<int>> posInPath(N);
	for (int s : order)
	{
		while (!dA[s].empty())
		{
			P.clear(); T.clear();
			P.push_back(s); posInPath[s].push_back(0);
			bool needA = true;
			while (true)
			{
				int v = P.back();
				std::vector<int> & lst = needA ? dA[v] : dB[v];
				if (lst.empty()) break; // stuck (degree mismatch at the depot): the path is discarded
				int idx = (int)(params.ran() % lst.size());
				int w = lst[idx];
				lst[idx] = lst.back(); lst.pop_back();
				removeOne(needA ? dA[w] : dB[w], v);
				P.push_back(w); T.push_back(needA ? 1 : 0);
				// closure: an earlier occurrence of w whose outgoing path edge has the other type
				int found = -1;
				for (int k = (int)posInPath[w].size() - 1; k >= 0; k--)
				{
					int i = posInPath[w][k];
					if (T[i] != T.back()) { found = i; break; }
				}
				posInPath[w].push_back((int)P.size() - 1);
				if (found >= 0)
				{
					std::vector<ABEdge> cyc;
					for (int i = found; i < (int)T.size(); i++) cyc.push_back({ P[i], P[i + 1], T[i] == 1 });
					cycles.push_back(cyc);
					for (int i = found + 1; i < (int)P.size(); i++) posInPath[P[i]].pop_back();
					P.resize(found + 1);
					T.resize(found);
					if (P.size() == 1) break;
					needA = (T.back() == 0);
				}
				else needA = !needA;
			}
			for (int x : P) posInPath[x].clear();
		}
	}
	if (cycles.empty()) return false;

	// Localized block of AB-cycles: a random seed cycle plus cycles touching the neighbourhood of the selection
	std::vector<char> selected(cycles.size(), 0);
	std::vector<char> nearSel(N, 0);
	std::vector<std::vector<int>> cyclesOfNode(N);
	for (int c = 0; c < (int)cycles.size(); c++)
		for (const ABEdge & e : cycles[c]) if (e.u != 0) cyclesOfNode[e.u].push_back(c);
	int totalEdges = 0;
	for (const auto & c : cycles) totalEdges += (int)c.size();
	int target = 4 + (int)(params.ran() % std::max<int>(1, std::min<int>(totalEdges / 2, 30)));
	int c0 = (int)(params.ran() % cycles.size());
	std::vector<int> frontier;
	int selEdges = 0;
	auto selectCycle = [&](int c)
	{
		selected[c] = 1;
		selEdges += (int)cycles[c].size();
		for (const ABEdge & e : cycles[c])
			if (e.u != 0 && !nearSel[e.u])
			{
				nearSel[e.u] = 1;
				for (int x : cyclesOfNode[e.u]) if (!selected[x]) frontier.push_back(x);
				for (int k = 0; k < std::min<int>(5, (int)params.correlatedVertices[e.u].size()); k++)
				{
					int y = params.correlatedVertices[e.u][k];
					for (int x : cyclesOfNode[y]) if (!selected[x]) frontier.push_back(x);
				}
			}
	};
	selectCycle(c0);
	while (selEdges < target && !frontier.empty())
	{
		int idx = (int)(params.ran() % frontier.size());
		int c = frontier[idx];
		frontier[idx] = frontier.back(); frontier.pop_back();
		if (!selected[c]) selectCycle(c);
	}

	// Offspring edge set: A minus selected A-edges plus selected B-edges
	std::vector<std::vector<int>> adj = adjA;
	for (int c = 0; c < (int)cycles.size(); c++)
		if (selected[c])
			for (const ABEdge & e : cycles[c])
			{
				if (e.isA) { removeOne(adj[e.u], e.v); removeOne(adj[e.v], e.u); }
				else { adj[e.u].push_back(e.v); adj[e.v].push_back(e.u); }
			}

	// Extract the routes (walks from the depot) and the subtours
	std::vector<std::vector<int>> routes, subtours;
	std::vector<char> visited(N, 0);
	while (!adj[0].empty())
	{
		int w = adj[0].back(); adj[0].pop_back();
		removeOne(adj[w], 0);
		std::vector<int> route;
		int prev = 0, cur = w;
		while (cur != 0)
		{
			if (visited[cur] || adj[cur].size() != 1) return false; // safety
			visited[cur] = 1;
			route.push_back(cur);
			int nxt = adj[cur][0];
			adj[cur].clear();
			removeOne(adj[nxt], cur);
			prev = cur; cur = nxt;
		}
		(void)prev;
		routes.push_back(route);
	}
	for (int v = 1; v < N; v++)
	{
		if (visited[v]) continue;
		if (adj[v].size() != 2) return false; // safety
		std::vector<int> st;
		int cur = v;
		while (!visited[cur])
		{
			visited[cur] = 1;
			st.push_back(cur);
			int nxt = -1;
			for (int x : adj[cur]) if (!visited[x]) { nxt = x; break; }
			if (nxt < 0) break;
			cur = nxt;
		}
		subtours.push_back(st);
	}

	// Merge each subtour into a route by the cheapest 2-opt-style reconnection with a nearby route edge
	std::vector<int> routeOf(N, -1), posOf(N, -1);
	auto indexRoute = [&](int r) { for (int p = 0; p < (int)routes[r].size(); p++) { routeOf[routes[r][p]] = r; posOf[routes[r][p]] = p; } };
	for (int r = 0; r < (int)routes.size(); r++) indexRoute(r);
	const auto & d = params.timeCost;
	std::vector<char> merged(subtours.size(), 0);
	int nbLeft = (int)subtours.size();
	for (int pass = 0; pass < 3 && nbLeft > 0; pass++)
	{
		for (int s = 0; s < (int)subtours.size(); s++)
		{
			if (merged[s]) continue;
			const std::vector<int> & st = subtours[s];
			int L = (int)st.size();
			double best = 1.e30; int bR = -1, bPos = -1, bQ = -1; bool bForward = true;
			for (int q = 0; q < L; q++)
			{
				int u = st[q], v = st[(q + 1) % L];
				for (int x : params.correlatedVertices[u])
				{
					int r = routeOf[x];
					if (r < 0) continue;
					int p = posOf[x];
					for (int side = 0; side < 2; side++)
					{
						// insertion between a=routes[r][pa-1] (or depot) and b=routes[r][pa] (or depot), pa in [0, size]
						int pa = (side == 0) ? p : p + 1;
						int a = (pa == 0) ? 0 : routes[r][pa - 1];
						int b = (pa == (int)routes[r].size()) ? 0 : routes[r][pa];
						double base = - d[a][b] - d[u][v];
						double c1 = base + d[a][u] + d[v][b]; // a, u, u-1, ..., v, b (backwards)
						double c2 = base + d[a][v] + d[u][b]; // a, v, v+1, ..., u, b (forwards)
						if (c1 < best) { best = c1; bR = r; bPos = pa; bQ = q; bForward = false; }
						if (c2 < best) { best = c2; bR = r; bPos = pa; bQ = q; bForward = true; }
					}
				}
			}
			if (bR < 0 && pass < 2) continue;
			std::vector<int> seq;
			if (bR < 0)
			{
				// No route nearby: open the subtour as a new route at its cheapest edge
				double bc = 1.e30; int bq = 0;
				for (int q = 0; q < L; q++)
				{
					int u = st[q], v = st[(q + 1) % L];
					double c = d[0][v] + d[u][0] - d[u][v];
					if (c < bc) { bc = c; bq = q; }
				}
				for (int k = 1; k <= L; k++) seq.push_back(st[(bq + k) % L]);
				routes.push_back(seq);
				indexRoute((int)routes.size() - 1);
			}
			else
			{
				if (bForward) for (int k = 1; k <= L; k++) seq.push_back(st[(bQ + k) % L]);   // v ... u
				else for (int k = 0; k < L; k++) seq.push_back(st[((bQ - k) % L + L) % L]);   // u ... v
				routes[bR].insert(routes[bR].begin() + bPos, seq.begin(), seq.end());
				indexRoute(bR);
			}
			merged[s] = 1; nbLeft--;
		}
	}
	if ((int)routes.size() > params.nbVehicles) return false;

	// Store the offspring
	int idx = 0;
	for (int r = 0; r < params.nbVehicles; r++)
	{
		if (r < (int)routes.size()) result.chromR[r] = routes[r];
		else result.chromR[r].clear();
		if (r < (int)routes.size()) for (int c : routes[r]) result.chromT[idx++] = c;
	}
	if (idx != params.nbClients) return false; // safety
	result.evaluateCompleteCost(params);
	return true;
}
