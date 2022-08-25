#pragma once
#include "SingleParticleMobility.h"


class ParticleSystem : public VectorFieldSum
{
public:

	std::vector<SingleParticle> particles;
	VectorFieldSum traction;

	ParticleSystem(std::vector<RectCoord> cs = std::vector<RectCoord>(), int n = 0) 
	{
		particles.resize(cs.size());
		for (int i = 0; i < particles.size(); i++)
		{
			particles[i] = SingleParticle(cs[i], n);

			append(&particles[i].flow);
			traction.append(&particles[i].traction);
		}
	}

	ParticleSystem(std::vector<SphereData> ds, std::vector<RectCoord> cs, int nst, int nn)
	{

		if (ds.size() != cs.size())
			std::cout << "Error: Data size and centers size in paticle system are not equal.\n";

		particles.resize(cs.size());
		for (int i = 0; i < particles.size(); i++)
		{
			particles[i] = SingleParticle(ds[i], cs[i], nst, nn);

			append(&particles[i].flow);
			traction.append(&particles[i].traction);
		}
	}


	ParticleSystem(const ParticleSystem& ps)
	{
		particles = ps.particles;

		for (SingleParticle p : particles)
		{
			append(&p.flow);
			traction.append(&p.traction);
		}
	}

	ParticleSystem(const int& n)
	{
		particles.resize(n);
		for (SingleParticle p : particles)
		{
			append(&p.flow);
			traction.append(&p.traction);
		}
	}

	void push_back(SingleParticle p)
	{
		particles.push_back(p);
		append(&particles.back().flow);
		traction.append(&particles.back().traction);
	}

	ParticleSystem operator +(const ParticleSystem& ps)
	{
		ParticleSystem temp(*this);

		if (particles.size() != ps.particles.size())
			std::cout << "Error:Summing different numbers of particles.\n";

		for (int i = 0; i < particles.size(); i++)
			temp.particles[i] = particles[i] + ps.particles[i];

		return temp;

	}

	ParticleSystem operator -(const ParticleSystem& ps)
	{
		ParticleSystem temp(*this);

		if (particles.size() != ps.particles.size())
			std::cout << "Error:Summing different numbers of particles.\n";

		for (int i = 0; i < particles.size(); i++)
			temp.particles[i] = particles[i] - ps.particles[i];

		return temp;

	}

	ParticleSystem operator +=(const ParticleSystem& ps)
	{
		*this = *this + ps;
		return *this;
	}

	ParticleSystem operator *=(const double& a)
	{
		*this = (*this) * a;
		return *this;
	}

	ParticleSystem operator *(const double& a)
	{
		ParticleSystem temp(*this);


		for (int i = 0; i < particles.size(); i++)
			temp.particles[i] = particles[i] *a;

		return temp;

	}

	ParticleSystem operator =(const ParticleSystem& ps)
	{
		this->particles = ps.particles;

		for (SingleParticle p : particles)
		{
			append(&p.flow);
			traction.append(&p.traction);
		}

		return *this;

	}

	double norm()
	{
		double temp = 0.0;
		for (SingleParticle p : particles)
		{
			double normj = p.norm();
			temp += normj * normj;
		}
		return sqrt(temp);

	}

	double dot(ParticleSystem& ps)
	{
		double temp = 0.0;
		for (int i =0; i < particles.size(); i++)
		{
			temp += particles[i].dot(ps.particles[i]);
		}
		return temp;

	}

	static double dot(ParticleSystem& particle, ParticleSystem& qarticle)
	{
		return particle.dot(qarticle);
	}

	int getN()
	{
		return particles.size();
	}

	void axpy(ParticleSystem* ps, double scale)
	{
		*this = *this + (*ps) * scale;
	}

	


};


class PSLHSOperator
{
public:
	ParticleSystem operator *(ParticleSystem& ps)
	{
		LHSOperator L;

		ParticleSystem temp(ps.getN());

		//persorm self interation part of operator.
		for (int i = 0; i < ps.particles.size(); i++)
			temp.particles[i] = L * ps.particles[i];

		for (int i = 0; i < ps.particles.size(); i++)
		{
			for (int j = 0; j < ps.particles.size(); j++)
			{
				if (i != j)
					for (int t = 0; t < NUMGLNODES; t++)
						for (int p = 0; p < temp.particles[i].numNodes; p++)
						{

							SurfaceCoord s(PI / 2.0 * (GLnodes[t] + 1), 2.0 * PI * (double)p / (double)ps.particles[i].numNodes);
							SphereCoord x(s, ps.particles[i].center);
							temp.particles[i].data[0][t][p] += ps.particles[j].traction(x);
						}

			}
			temp.particles[i].refreshData();
			std::cout << "Operator on particle " << i + 1 << " computed.\n";
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

			//persorm self interation part of operator.
			for (int i = 0; i < ps.particles.size(); i++)
				temp.particles[i] = R * ps.particles[i];

			for (int i = 0; i < ps.particles.size(); i++)
			{
				for (int j = 0; j < ps.particles.size(); j++)
				{
					if (i != j)
						for (int t = 0; t < NUMGLNODES; t++)
							for (int p = 0; p < temp.particles[j].numNodes; p++)
							{

								SurfaceCoord s(PI / 2.0 * (GLnodes[t] + 1), 2.0 * PI * (double)p / (double)ps.particles[i].numNodes);
								SphereCoord x(s, ps.particles[i].center);
								temp.particles[i].data[0][t][p] -= ps.particles[j].traction(x);
							}


				}
				temp.particles[i].refreshData();
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






ParticleSystem SolveMobilityManySphere(std::vector<RectCoord> cs,std::vector<RectCoord> Fs, std::vector<RectCoord> Ts, int seriesSize)
{

	SphericalVectorField initialguess(&e_r);

	if ((cs.size() != Fs.size()) || (Fs.size() != Ts.size()))
	{
		std::cout << "Error: size mismatch between mobility problem input vectors.\n";
		return ParticleSystem();
	}

	int size = cs.size();

	ParticleSystem BIEsolution;
	for(int i = 0; i < size; i++)
		BIEsolution.push_back(SingleParticle(discretize(&initialguess), cs[i], seriesSize, NUMTRAPNODES));

	PSLHSOperator L;

	std::vector<ForceBalance> rhos;

	for (int i = 0; i < size; i++)
		rhos.push_back(ForceBalance(Fs[i], Ts[i]));


	ParticleSystem rh;
	for (int i = 0; i < size; i++)
	{
		SphereData* temp = discretize(&rhos[i]);
		rh.push_back(SingleParticle(temp, cs[i], seriesSize, NUMTRAPNODES));
		delete temp;
	}

	PSRHSOperator R;
	ParticleSystem rhs = R * rh;

	PSIdentityPreconditioner* I = new PSIdentityPreconditioner();

	GMRES(&L, &BIEsolution, &rhs, I, 20, 5, 1e-6);



	ParticleSystem soln = BIEsolution + rh;

	delete I;
	return soln;

}