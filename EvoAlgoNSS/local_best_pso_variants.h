#pragma once

#include "local_best_pso.h"

template<typename T>
struct PSOstruct_inertia : PSOstruct<T>
{
	PSOstruct_inertia(const T& i_c1, const T& i_c2, const size_t& i_sneigh, const T& i_w, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: PSOstruct{ i_c1 , i_c2 , i_sneigh, i_npop, i_tol, i_iter_max }, w{ i_w }
	{
		assert(w > 0);
	}
	T w;
};

template<typename T>
struct PSOstruct_clamping : PSOstruct<T>
{
	PSOstruct_clamping(const T& i_c1, const T& i_c2, const size_t& i_sneigh, const T& i_alpha, const std::vector<T>& i_vmax, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: PSOstruct{ i_c1 , i_c2 , i_sneigh, i_npop, i_tol, i_iter_max }, alpha{ i_alpha }, vmax{ i_vmax }
	{
		assert(alpha > 0);
		for (const auto& p : vmax) { assert(p > 0); };
	}
	T alpha;
	std::vector<T> vmax;
};

template<typename T, typename F>
class Solver<T, F, PSOstruct_inertia<T>> final : public Solver<T, F, PSOstruct<T>>
{
public:
	Solver(const PSOstruct_inertia<T>& pso, const Population<T>& popul) : Solver < T, F, PSOstruct<T>>{ { pso.c1, pso.c2, pso.sneigh, pso.npop, pso.tol, pso.iter_max}, popul }
	{
		w = pso.w;
	}
private:
	T w;
	void velocity_update(const size_t& iter);
};

template<typename T, typename F>
void Solver<T, F, PSOstruct_inertia<T>>::velocity_update(const size_t& iter)
{
	const auto& r = generate_r(ndv);
	for (auto k = 0; k < nneigh; ++k)
	{
		for (auto l = 0; l < neighbourhoods[k].size(); ++l)
		{
			auto i = neighbourhoods[k][l];
			for (auto j = 0; j < ndv; ++j)
			{
				{
					velocity[i][j] = w * velocity[i][j] + c1 * r[0][j] * (personal_best[i][j]
						- individuals[i][j]) + c2 * r[1][j] * (local_best[k][j] - individuals[i][j]);
					individuals[i][j] = individuals[i][j] + velocity[i][j];
				}
			}
		}
	}
	w = ((w - 0.4) * (iter_max - iter)) / (iter_max + 0.4);
}

template<typename T, typename F>
class Solver<T, F, PSOstruct_clamping<T>> final : public Solver<T, F, PSOstruct<T>>
{
public:
	Solver(const PSOstruct_clamping<T>& pso, const Population<T>& popul) : Solver < T, F, PSOstruct<T>>{ { pso.c1, pso.c2, pso.sneigh, pso.npop, pso.tol, pso.iter_max}, popul }
	{
		alpha = pso.alpha;
		vmax = pso.vmax;
		assert(vmax.size() == ndv);
	}
private:
	T alpha;
	std::vector<T> vmax;
	void velocity_update(const size_t& iter);
};

template<typename T, typename F>
void Solver<T, F, PSOstruct_clamping<T>>::velocity_update(const size_t& iter)
{
	const auto& r = generate_r(ndv);
	for (auto k = 0; k < nneigh; ++k)
	{
		for (auto l = 0; l < neighbourhoods[k].size(); ++l)
		{
			auto i = neighbourhoods[k][l];
			for (auto j = 0; j < ndv; ++j)
			{
				{
					velocity[i][j] = velocity[i][j] + c1 * r[0][j] * (personal_best[i][j]
						- individuals[i][j]) + c2 * r[1][j] * (local_best[k][j] - individuals[i][j]);
					if (velocity[i][j] > vmax[j])
					{
						velocity[i][j] = vmax[j];
					}
					individuals[i][j] = individuals[i][j] + velocity[i][j];
				}
			}
		}
	}
	for (auto& p : vmax)
	{
		p = (1 - std::pow(iter / iter_max, alpha)) * p;
	}
}