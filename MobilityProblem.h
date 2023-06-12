#pragma once
#include "SingleParticleMobility.h"
#include <ctime>

class ParticleSystem : public SphericalVectorField
{
public:

	std::vector<SingleParticle> particles;
	VectorFieldSum traction;

	ParticleSystem(std::vector<RectCoord> cs = std::vector<RectCoord>(), int n = 6) 
	{
		particles.resize(cs.size());
		for (size_t i = 0; i < particles.size(); i++)
		{
			particles[i] = SingleParticle(cs[i], n);

			traction.append(&particles[i].traction);
		}
	}

	ParticleSystem(std::vector<SphereData> ds, std::vector<RectCoord> cs, int nst)
	{

		if (ds.size() != cs.size())
			std::cout << "Error: Data size and centers size in paticle system are not equal.\n";

		particles.resize(cs.size());
		for (size_t i = 0; i < particles.size(); i++)
		{
			particles[i] = SingleParticle(ds[i], cs[i], nst);

			traction.append(&particles[i].traction);
		}

	
	}


	ParticleSystem(const ParticleSystem& ps)
	{
		particles = ps.particles;

		traction.clear();

		for (int i = 0; i < particles.size(); i++)
		{
			traction.append(&particles[i].traction);
		}
	}


	ParticleSystem(const int& n)
	{
		particles.resize(n);
		for (int i = 0; i < particles.size(); i++)
		{
			traction.append(&particles[i].traction);
		}
	}

	void push_back(const SingleParticle & p)
	{
		particles.push_back(p);
		traction.clear();
		for (int i = 0; i < particles.size(); i++)
		{
			traction.append(&particles[i].traction);
		}
	}

	ParticleSystem operator +(const ParticleSystem& ps) const
	{
		ParticleSystem temp(*this);

		if (particles.size() != ps.particles.size())
			std::cout << "Error:Summing different numbers of particles.\n";

		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] = particles[i] + ps.particles[i];

		return temp;

	}

	ParticleSystem operator -(const ParticleSystem& ps) const
	{
		ParticleSystem temp(*this);

		if (particles.size() != ps.particles.size())
			std::cout << "Error:Summing different numbers of particles.\n";

		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] = particles[i] - ps.particles[i];

		return temp;

	}

	ParticleSystem operator -() const
	{
		ParticleSystem temp(*this);


		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] = -particles[i];

		return temp;

	}

	ParticleSystem & operator +=(const ParticleSystem& ps)
	{
		for (int i = 0; i < particles.size(); i++)
			particles[i] += ps.particles[i];

		return *this;
	}

	ParticleSystem & operator *=(const double& a)
	{
		for (int i = 0; i < particles.size(); i++)
			particles[i] *= a;
		
		return *this;
	}

	ParticleSystem operator *(const double& a) const
	{
		ParticleSystem temp(*this);


		for (size_t i = 0; i < particles.size(); i++)
			temp.particles[i] *=a ;

		return temp;

	}

	ParticleSystem& operator =(const ParticleSystem& ps)
	{
		particles = ps.particles;

		traction.clear();

		for (int i = 0; i < particles.size(); i++)
		{
				traction.append(&particles[i].traction);
		}

		return *this;

	}

	


	double norm() const
	{
		double temp = 0.0;
		for (int i = 0; i < particles.size();i++)
		{
			double normj = particles[i].norm();
			temp += normj * normj;
		}
		return sqrt(temp);

	}

	double dot(const ParticleSystem& ps) const
	{
		double temp = 0.0;
		for (size_t i =0; i < particles.size(); i++)
		{
			temp += particles[i].dot(ps.particles[i]);
		}
		return temp;

	}

	static double dot(const ParticleSystem& particles, const ParticleSystem& qarticles)
	{
		return particles.dot(qarticles);
	}

	int getN() const
	{
		return particles.size();
	}

	void axpy(ParticleSystem* ps, double scale)
	{
		for (int i = 0; i < particles.size(); i++)
			particles[i].axpy(&(ps->particles[i]), scale);
		
	}


	RectCoord operator()(const SphereCoord& x) const
	{
		RectCoord temp;

		for (int i = 0; i < particles.size(); i++)
		{
			RectCoord r = x.center - particles[i].center;
			if ( sqrt(r.x * r.x + r.y * r.y + r.z * r.z) < 4.0)
				temp += particles[i].flow(x);
			else
				temp += particles[i].quadflow(x);
		}

		return temp;
	}

};


class PSLHSOperator
{
public:
	ParticleSystem operator *(const ParticleSystem& ps)
	{
		LHSOperator L;

		ParticleSystem temp(ps.getN());

		for (size_t i = 0; i < ps.particles.size(); i++)
		{
			std::time_t ti = std::time(0);
			temp.particles[i] = L * ps.particles[i];
			RectCoord center = ps.particles[i].center;
			SingleParticle tr(center, ps.particles[i].numSeriesTerms);
			
			SphereData t(2 * ps.particles[i].numSeriesTerms + 1 , 2* ps.particles[i].numSeriesTerms + 1);

			

			for (size_t j = 0; j < ps.particles.size(); j++)
			{
				if (i != j)
				{
					if(norm(ps.particles[i].center - ps.particles[j].center) < 4.0)
						t = t + discretize(&ps.particles[j].traction, center , ps.particles[j].consts);
					else
						t = t + discretize(&ps.particles[j].quadtraction, center , ps.particles[j].consts);
				}
			}
			temp.particles[i] += SingleParticle(t , center , ps.particles[i].numSeriesTerms);
			std::cout << "Operator on particle " << i + 1 << " computed. Time elapsed: " << std::time(0) - ti <<"\n";
		}
		return temp;
	}
};

	class PSRHSOperator
	{
	public:
		ParticleSystem operator *(ParticleSystem& ps)
		{
			RHSOperator R;

			ParticleSystem temp(ps.getN());

			for (size_t i = 0; i < ps.particles.size(); i++)
			{
				temp.particles[i] = R * ps.particles[i];
				RectCoord center = ps.particles[i].center;
				SingleParticle tr(center, ps.particles[i].numSeriesTerms);


				SphereData t(2 * ps.particles[i].numSeriesTerms + 1, 2 * ps.particles[i].numSeriesTerms + 1);

				for (size_t j = 0; j < ps.particles.size(); j++)
				{
					if (i != j)
					{
						if (norm(ps.particles[i].center - ps.particles[j].center) < 4.0)
							t = t - discretize(&ps.particles[j].traction, center , ps.particles[j].consts);
						else
							t = t - discretize(&ps.particles[j].quadtraction, center, ps.particles[j].consts);
					}
				}
				temp.particles[i] += SingleParticle(t, center, ps.particles[i].numSeriesTerms);
				std::cout << "RHS Operator on particle " << i + 1 << " computed.\n";
			
			}
			return temp;
		}

		
	};
class PSIdentityPreconditioner
{
public:
	ParticleSystem solve(const ParticleSystem& ps)
	{
		return ps;
	}
};


//Solve single particle problem with linear combinations of F,T, compare to exact solution off sphere.
// evaulate solution for single sphere problem off sphere, compare to exact solution off sphere.
//Two particles widely seperated, different forces and torques, compare solutions to isolated sphere problem.
// Show difference in data using sup norm.

ParticleSystem SolveMobilityManySphere(std::vector<RectCoord> cs,std::vector<RectCoord> Fs, std::vector<RectCoord> Ts, int seriesSize ,const ParticleSystem & initial= ParticleSystem(0), int numIters = 100, double r = 1e-9)
{

	if ((cs.size() != Fs.size()) || (Fs.size() != Ts.size()))
	{
		std::cout << "Error: size mismatch between mobility problem input vectors.\n";
		return ParticleSystem(0);
	}

	int size = cs.size();
	ParticleSystem BIEsolution(cs , seriesSize);
	if (initial.particles.size() == size)
		BIEsolution = initial;
	
		
	PSLHSOperator L;



    ParticleSystem rh(cs , seriesSize);
	
	for (int i = 0; i < size; i++)
	{
		ForceBalance rho = ForceBalance(Fs[i], Ts[i]);
		rh.particles[i] = SingleParticle(discretize(&rho , cs[i] , rh.particles[i].consts), cs[i], seriesSize);
	}


	BIEsolution;
	PSRHSOperator R;
	ParticleSystem rhs = R * rh;

	PSIdentityPreconditioner I;

	GMRES(&L, &BIEsolution, &rhs, &I, numIters, 1, r);



	ParticleSystem soln = BIEsolution + rh;

	//std::cout << Integrate(&BIEsolution, 1.0, BIEsolution.particles[0].center) * 0.25 / PI << std::endl;
	//std::cout << Integrate(&rh, 1.0, rh.particles[0].center) * 0.25 / PI << std::endl;

	return soln;

}