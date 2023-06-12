#pragma once
#include "StokesOperators.h"
#include "GMRES.h"

const int MAXITER = 100;
int numcount = 0;
class SingleParticle
{
public:
	int numSeriesTerms;

	RectCoord center;

	
	SphereData data;
	VSHSeries fourierdata;


	MathConstants consts;
	SingleLayerDirectDiscrete quadflow;
	StokesTractionDirectDiscrete quadtraction;

	SingleLayerHarmonic flow;
	StokesPressure pressure;
	StokesTractionHarmonic traction;
	

	SingleParticle(RectCoord c = RectCoord(),int n = 6  ) :
		  data(2*n + 1,2*n + 1),
		  consts(2 * n+1 , 2 * n + 1),
		  fourierdata(n, &consts, c),
		  quadflow(&data,&consts, c),
		  quadtraction(&data,&consts , c),
		  flow(n , c) ,
		  pressure(n,c) , 
		  traction(&flow , &pressure),
		  center(c),
		  numSeriesTerms(n)
	{
		fourierdata.approximate(data, numSeriesTerms);
		flow.solve(fourierdata);
		pressure.solve(fourierdata);
	}

	SingleParticle(const SphereData& d, RectCoord c, int nst) :
		data(d),
		consts(2 * nst + 1, 2 * nst + 1),
		numSeriesTerms(nst),
		fourierdata(nst, &consts, c),
		quadflow(&data, &consts, c),
		quadtraction(&data, &consts, c),
		flow(nst , c),
		pressure(nst , c),
		traction(&flow, &pressure),
		center(c)
	{
		if (data.size() != 2 * nst + 1)
			std::cout << "Error! size mismatch in dimensions of particle.\n";
		fourierdata.approximate(data, numSeriesTerms);
		flow.solve(fourierdata);
		pressure.solve(fourierdata);
	}


	

	SingleParticle(const SingleParticle& particle) :
		numSeriesTerms(particle.numSeriesTerms) ,
		data(particle.data),
		consts(2*particle.numSeriesTerms + 1 , 2 * particle.numSeriesTerms + 1),
		quadflow(&data,&consts, particle.center),
		quadtraction(&data,&consts, particle.center),
		fourierdata(particle.fourierdata),
		flow(particle.flow), 
		pressure(particle.pressure),
		traction(&flow, &pressure),
		center(particle.center)
	{
	}

	SingleParticle operator +(const SingleParticle & particle) const
	{
		SingleParticle temp(center , numSeriesTerms);

		
			for (int t = 0; t < (int)consts.NUMGLNODES; t++)
				for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
					temp.data[t][p] = data[t][p] + particle.data[t][p];
		

		temp.fourierdata = fourierdata + particle.fourierdata;


		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.center = center;

		return temp;

	}

    SingleParticle operator -(const SingleParticle & particle) const
	{
		SingleParticle temp(center , numSeriesTerms);

		for (int t = 0; t < consts.NUMGLNODES; t++)
			for (int p = 0; p < consts.NUMTRAPNODES; p++)
				temp.data[t][p] = data[t][p] - particle.data[t][p];

		temp.fourierdata = fourierdata - particle.fourierdata;

		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.center = center;
		return temp;

	}

	SingleParticle operator -() const
	{
		SingleParticle temp(center , numSeriesTerms);

		for (int t = 0; t < consts.NUMGLNODES; t++)
			for (int p = 0; p < consts.NUMTRAPNODES; p++)
				temp.data[t][p] = - data[t][p];

		temp.fourierdata = -fourierdata;

		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(temp.fourierdata);
		temp.pressure.solve(temp.fourierdata);

		temp.center = center;
		return temp;
	}

    SingleParticle & operator +=(const SingleParticle & particle)
	{
		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				data[t][p] += particle.data[t][p];


		fourierdata += particle.fourierdata;



		flow.solve(fourierdata);
		pressure.solve(fourierdata);
		return *this;
	}

	SingleParticle & operator *=(const double & a)
	{
		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				data[t][p] *= a;


		fourierdata *= a;



		flow.solve(fourierdata);
		pressure.solve(fourierdata);
		return *this;
	}

	SingleParticle operator *(const double& a) const
	{
		SingleParticle temp(center , numSeriesTerms);
		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				temp.data[t][p] = data[t][p] * a;

		temp.fourierdata = fourierdata * a;

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
		consts = particle.consts;
		flow = particle.flow;
		pressure = particle.pressure;
		center = particle.center;

		return *this;

	}



	double norm() const
	{
		return sqrt(dot(*this));
	}

	double dot(const SingleParticle& particle) const
	{
		return L2InnerProductDiscrete(data , particle.data );
	}

	static double dot(const SingleParticle &particle, const SingleParticle& qarticle)
	{
		return particle.dot(qarticle);
	}

	int getN() const
	{
		return (int)consts.NUMGLNODES * (int)consts.NUMTRAPNODES * 3;
	}

	void axpy(SingleParticle* particle, double scale)
	{


		for (int t = 0; t < (int)consts.NUMGLNODES; t++)
			for (int p = 0; p < (int)consts.NUMTRAPNODES; p++)
				data[t][p] +=particle->data[t][p] * scale;


		fourierdata.axpy( particle->fourierdata , scale);



		flow.solve(fourierdata);
		pressure.solve(fourierdata);

	}

	

	void refreshData()
	{
		fourierdata.approximate(data, numSeriesTerms);
		

		flow.solve(fourierdata);
		pressure.solve(fourierdata);
	}

	
};

RectCoord rotationCoefficient(const SphereData& data , const RectCoord & center,const double & radius , const MathConstants & consts)
{
	

		RectCoord total = 0.0;

		MathConstants temp;
		if (consts.NUMGLNODES != data.size() || consts.NUMTRAPNODES != data[0].size())
			temp = MathConstants(data[0].size(), data.size());
		else
			temp = consts;

		for (int p = 0; p < temp.NUMTRAPNODES; p++)
			for (int i = 0; i < temp.NUMGLNODES; i++)
			{
				SurfaceCoord s(temp.PI / 2.0 * (temp.GLnodes[i] + 1), 2.0 * temp.PI * (double)p / (double)temp.NUMTRAPNODES);
				SphereCoord x(radius, s, center);
				total = total + radius * radius * sin(s.theta) * temp.GLweights[i] * cross(SphereToRect(x) - x.center,data[i][p]) * (temp.PI / temp.NUMTRAPNODES) * consts.PI;


			}

		return total;
	
}

class LHSOperator
{
public:
	SingleParticle operator *(const SingleParticle & particle) const
	{
		SingleParticle temp(particle.center , particle.numSeriesTerms);

		RectCoord average = 0.25 / particle.consts.PI * IntegrateDiscrete(particle.data , 1.0,particle.center, particle.consts );
		RectCoord rcoef = 0.375/particle.consts.PI * rotationCoefficient(particle.data, particle.center, 1.0 , particle.consts);

		for(int t = 0; t < particle.consts.NUMGLNODES; t++)
			for (int p = 0; p < particle.consts.NUMTRAPNODES; p++)
			{
				SurfaceCoord s(particle.consts.PI / 2.0 * (particle.consts.GLnodes[t] + 1), 2.0 * particle.consts.PI * (double)p / (double)particle.consts.NUMTRAPNODES);
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


		for (int t = 0; t < particle.consts.NUMGLNODES; t++)
		{
			for (int p = 0; p < particle.consts.NUMTRAPNODES; p++)
			{
				SurfaceCoord s(particle.consts.PI / 2.0 * (particle.consts.GLnodes[t] + 1), 2.0 * particle.consts.PI * (double)p / (double)particle.consts.NUMTRAPNODES);
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
	double PI = 4.0 * atan(1.0);
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
	
	const double PI = 4.0 * atan(1.0);

	MathConstants consts(2 * seriesSize + 1, 2 * seriesSize + 1);

	SingleParticle rh(discretize(&rho, RectCoord() , consts ), RectCoord(), seriesSize);
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
	std::cout << "Angular Velocity of rho: " << rotationCoefficient(discretize(&rh.flow), RectCoord(), 1.0 , rh.consts) * 0.375 / PI << "\n";
	return soln;

}

SingleParticle SolveBIEWithKnown(int seriesSize)
{



	LHSOperator L;

	std::default_random_engine gen;
	std::uniform_real_distribution rand;

	SphereData data(2 * seriesSize + 1, 2 * seriesSize + 1);

	

	for (size_t t = 0; t < data.size(); t++)
	{
		for (size_t p = 0; p < data[t].size(); p++)
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
	for (int t = 0; t < err.consts.NUMGLNODES - 1; t++)
	{
		std::cout << "[";
		for (int p = 0; p < err.consts.NUMTRAPNODES - 1; p++)
		  std::cout << norm(err.data[t][p]) << ", ";
		std::cout << norm(err.data[t][err.consts.NUMTRAPNODES - 1]) << " ], ";
	}
	std::cout << "[";
	for (int p = 0; p < err.consts.NUMTRAPNODES - 1; p++)
		std::cout << norm(err.data[(int)(err.consts.NUMTRAPNODES - 1)][p]) << ", ";
	std::cout << norm(err.data[(int)(err.consts.NUMGLNODES - 1)][err.consts.NUMTRAPNODES - 1]) <<" ]]" << std::endl;


	return solution;

}