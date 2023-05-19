#pragma once
#include "StokesOperators.h"
#include "GMRES.h"

const int MAXITER = 100;
class SingleParticle
{
public:
	int numSeriesTerms;

	RectCoord center;


	SphereData data;
	VSHSeries fourierdata;

	SingleLayerDirectDiscrete quadflow;
	StokesTractionDirectDiscrete quadtraction;

	SingleLayerHarmonic flow;
	StokesPressure pressure;
	StokesTractionHarmonic traction;
	

	SingleParticle(RectCoord c = RectCoord(),int n = 5  ) :
		  fourierdata(n , c),
		  quadflow(&data , c),
		  quadtraction(&data , c),
		  flow(n , c) ,
		  pressure(n,c) , 
		  traction(&flow , &pressure),
		  center(c),
		  numSeriesTerms(n)
	{
		data.resize((int)NUMGLNODES);
		for (int i = 0; i < (int)NUMGLNODES; i++)
		{
			data[i].resize((int)NUMTRAPNODES);
			for (int j = 0; j < NUMTRAPNODES; j++)
				data[i][j] = 0.0;
		}
		fourierdata.approximate(data, numSeriesTerms);
		flow.solve(fourierdata);
		pressure.solve(fourierdata);
	}

	SingleParticle(const SphereData& d, RectCoord c ,int nst):  
		numSeriesTerms(nst),
		data(d),
		quadflow(&data, c),
		quadtraction(&data, c),
		fourierdata(nst , c),
		flow(nst , c),
		pressure(nst , c),
		traction(&flow, &pressure),
		center(c)
	{
		fourierdata.approximate(data, numSeriesTerms);
		flow.solve(fourierdata);
		pressure.solve(fourierdata);
	}


	

	SingleParticle(const SingleParticle& particle) :
		numSeriesTerms(particle.numSeriesTerms) ,
		data(particle.data),
		quadflow(&data, particle.center),
		quadtraction(&data, particle.center),
		fourierdata(particle.fourierdata),
		flow(particle.flow), 
		pressure(particle.pressure),
		traction(&flow, &pressure),
		center(particle.center)
	{}

	SingleParticle(SingleParticle&& particle) noexcept:
		numSeriesTerms(particle.numSeriesTerms),
		data(std::move(particle.data)),
		quadflow(&data, particle.center),
		quadtraction(&data, particle.center),
		fourierdata(std::move(particle.fourierdata)),
		flow(std::move(particle.flow)),
		pressure(std::move(particle.pressure)),
		traction(&flow, &pressure),
		center(particle.center)
	{}

	SingleParticle operator +(const SingleParticle & particle)
	{
		SingleParticle temp(center , numSeriesTerms);

		
			for (int t = 0; t < (int)NUMGLNODES; t++)
				for (int p = 0; p < (int)NUMTRAPNODES; p++)
					temp.data[t][p] = data[t][p] + particle.data[t][p];
		

		temp.fourierdata = std::move(fourierdata + particle.fourierdata);


		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.center = center;

		return temp;

	}

    SingleParticle operator -(const SingleParticle & particle)
	{
		SingleParticle temp(center , numSeriesTerms);

		for (int t = 0; t < NUMGLNODES; t++)
			for (int p = 0; p < NUMTRAPNODES; p++)
				temp.data[t][p] = data[t][p] - particle.data[t][p];

		temp.fourierdata = std::move(fourierdata - particle.fourierdata);

		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.center = center;
		return temp;

	}

	SingleParticle operator -() const
	{
		SingleParticle temp(center , numSeriesTerms);

		for (int t = 0; t < NUMGLNODES; t++)
			for (int p = 0; p < NUMTRAPNODES; p++)
				temp.data[t][p] = - data[t][p];

		temp.fourierdata = std::move(-fourierdata);

		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.center = center;
		return temp;
	}

    SingleParticle & operator +=(const SingleParticle & particle)
	{
		*this = std::move(*this + particle);
		return *this;
	}

	SingleParticle & operator *=(const double & a)
	{
		*this = std::move((*this) * a);
		return *this;
	}

	SingleParticle operator *(const double& a)
	{
		SingleParticle temp(center , numSeriesTerms);
		for (int t = 0; t < (int)NUMGLNODES; t++)
			for (int p = 0; p < (int)NUMTRAPNODES; p++)
				temp.data[t][p] = data[t][p] * a;

		temp.fourierdata = std::move(fourierdata * a);

		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.center = center;

		return temp;
	}

	SingleParticle& operator =(const SingleParticle & particle)
	{
		data = particle.data;
		fourierdata = particle.fourierdata;
		numSeriesTerms = particle.numSeriesTerms;
		flow = particle.flow;
		pressure = particle.pressure;
		center = particle.center;

		return *this;

	}

	SingleParticle& operator =(SingleParticle&& particle) noexcept
	{
		data = std::move(particle.data);
		fourierdata = std::move(particle.fourierdata);
		numSeriesTerms = particle.numSeriesTerms;
		flow = std::move(particle.flow);
		pressure = std::move(particle.pressure);
		center = particle.center;

		return *this;

	}


	double norm()
	{
		return sqrt(dot(*this));
	}

	double dot(SingleParticle& particle)
	{
		return L2InnerProductDiscrete(data , particle.data );
	}

	static double dot(SingleParticle &particle, SingleParticle& qarticle)
	{
		return particle.dot(qarticle);
	}

	int getN()
	{
		return (int)NUMGLNODES * (int)NUMTRAPNODES * 3;
	}

	void axpy(SingleParticle* particle, double scale)
	{
		*this = *this + (*particle) * scale;
	}

	

	void refreshData()
	{
		VSHSeries newSeries(numSeriesTerms, center);
        newSeries.approximate(data , numSeriesTerms);
		fourierdata = newSeries;
		

		flow.solve(fourierdata);
		pressure.solve(fourierdata);
	}

	void operator ~()
	{

	}
};

RectCoord rotationCoefficient(const SphereData& data , const RectCoord & center,const double & radius)
{
	

		RectCoord total = 0.0;

		for (int p = 0; p < NUMTRAPNODES; p++)
			for (int i = 0; i < NUMGLNODES; i++)
			{
				SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
				SphereCoord x(radius, s, center);
				total = total + radius * radius * sin(s.theta) * GLweights[i] * cross(SphereToRect(x) - x.center,data[i][p]) * (PI / NUMTRAPNODES) * PI;


			}

		return total;
	
}

class LHSOperator
{
public:
	SingleParticle operator *( SingleParticle & particle)
	{
		SingleParticle temp(particle.center , particle.numSeriesTerms);

		RectCoord average = 0.25 / PI * IntegrateDiscrete(particle.data );
		RectCoord rcoef = 0.375/PI * rotationCoefficient(particle.data, particle.center, 1.0);

		for(int t = 0; t < NUMGLNODES; t++)
			for (int p = 0; p < NUMTRAPNODES; p++)
			{
				SurfaceCoord s(PI / 2.0 * (GLnodes[t] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
				SphereCoord x(s , particle.center);

				temp.data[t][p] =  particle.data[t][p] + particle.traction(x) + average + cross(rcoef, SphereToRect(x) - x.center);

			}

		temp.refreshData();

		return temp;
	}


};

class RHSOperator
{
public:
	SingleParticle operator *(const SingleParticle & particle)
	{
	    SingleParticle temp(particle.center , particle.numSeriesTerms);


		for (int t = 0; t < NUMGLNODES; t++)
		{
			for (int p = 0; p < NUMTRAPNODES; p++)
			{
				SurfaceCoord s(PI / 2.0 * (GLnodes[t] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
				SphereCoord x(s , particle.center);

				temp.data[t][p] =  -particle.data[t][p] - particle.traction(x);

			}
		}

		temp.refreshData();

		return temp;
	}
};

class IdentityPreconditioner
{
public:
	SingleParticle solve(const SingleParticle& particle)
	{
		return particle;
	}
};


class ForceBalance : public SphericalVectorField
{
private:
	RectCoord force;
	RectCoord torque;
public:

	ForceBalance(RectCoord f , RectCoord t) : force(f) , torque(t){}

	RectCoord operator()(const SphereCoord & x) const
	{
		return 0.25 / PI *  force + 0.375 / PI * cross(torque , SphereToRect(x) - x.center);
	}
};



SingleParticle SolveMobilitySingleSphere(RectCoord F , RectCoord T , int seriesSize)
{

		
	
	LHSOperator L;
	ForceBalance rho(F,T);
	

	SingleParticle rh(discretize(&rho), RectCoord(), seriesSize);
	RHSOperator R;
	SingleParticle rhs = R * rh;
	
	IdentityPreconditioner I;
    SingleParticle solution = -rh;
	
	GMRES(&L, &solution, &rhs, &I, 100, 5, 1e-10);


	
	SingleParticle soln = solution + rh;

	std::cout <<"Average of BIE solution: " << IntegrateDiscrete(solution.data) * 0.25 / PI << std::endl;
	std::cout <<"Average of rhs: " << IntegrateDiscrete(rh.data) * 0.25 / PI << std::endl;
	std::cout <<"Integral of flow: " <<Integrate(&soln.flow) * 0.25 / PI << std::endl;
	std::cout << "Integral of rho: " << Integrate(&rh.flow) * 0.25 / PI << std::endl;
	std::cout << "Angular Velocity of rho: " << rotationCoefficient(discretize(&rh.flow), RectCoord(), 1.0) * 0.375 / PI << "\n";
	return soln;

}

SingleParticle SolveBIEWithKnown(int seriesSize)
{



	LHSOperator L;

	std::default_random_engine gen;
	std::uniform_real_distribution rand;

	SphereData data;

	data.resize(NUMGLNODES);

	for (size_t t = 0; t < NUMGLNODES; t++)
	{
		data[t].resize(NUMTRAPNODES);
		for (size_t p = 0; p < NUMTRAPNODES; p++)
			data[t][p] = RectCoord(rand(gen), rand(gen), rand(gen));
	}
	SingleParticle rh(data, RectCoord(), seriesSize);
	
	SingleParticle rhs = L * rh;

	IdentityPreconditioner I;
	SingleParticle solution;

	GMRES(&L, &solution, &rhs, &I, 100, 5, 1e-10);

	SingleParticle err = solution - rh;

	std::cout << "Error in solving: " << err.norm()<< std::endl;

	
	std::cout << "[ ";
	for (int t = 0; t < NUMGLNODES - 1; t++)
	{
		std::cout << "[";
		for (int p = 0; p < NUMTRAPNODES - 1; p++)
		  std::cout << norm(err.data[t][p]) << ", ";
		std::cout << norm(err.data[t][NUMTRAPNODES - 1]) << " ], ";
	}
	std::cout << "[";
	for (int p = 0; p < NUMTRAPNODES - 1; p++)
		std::cout << norm(err.data[NUMTRAPNODES - 1][p]) << ", ";
	std::cout << norm(err.data[NUMGLNODES - 1][NUMTRAPNODES - 1]) <<" ]]" << std::endl;


	return solution;

}