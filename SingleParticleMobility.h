#pragma once
#include "StokesOperators.h"
#include "GMRES.h"

const int MAXITER = 100;

class SingleParticle
{
public:
	int numSeriesTerms;
	int numNodes;

	RectCoord center;


	SphereData* data;
	VSHSeries* fourierdata;
	SingleLayerHarmonic flow;
	StokesPressure pressure;
	StokesTractionHarmonic traction;

	SingleParticle(RectCoord c = RectCoord(),int n = 0  ) :  flow(5) , pressure(5) , traction(&flow , &pressure)
	{
		data = new SphereData();
		fourierdata = new VSHSeries(5, NUMTRAPNODES , c);


		for (auto i : *data)
			for (auto j : i)
				j = 0.0;

		numSeriesTerms = 5;
		numNodes = NUMTRAPNODES;

		(*fourierdata).approximate(*data, numSeriesTerms, numNodes);
		flow.solve(*fourierdata);
		pressure.solve(*fourierdata);
	}

	SingleParticle(const SphereData& d, RectCoord c ,int nst, int nn): numNodes(nn), numSeriesTerms(nst), flow(5), pressure(5), traction(&flow, &pressure)
	{
		data = new SphereData(d);
		fourierdata = new VSHSeries(nst, nn , c);

		fourierdata->approximate(*data, numSeriesTerms, numNodes);
		flow.solve(*fourierdata);
		pressure.solve(*fourierdata);
	}

	SingleParticle(SphereData* d, RectCoord c, int nst, int nn) : numNodes(nn), numSeriesTerms(nst), flow(nst), pressure(nst), traction(&flow, &pressure)
	{
		data = d;
		fourierdata = new VSHSeries(5, NUMTRAPNODES , c);

		(*fourierdata).approximate(*data, numSeriesTerms, numNodes);
		flow.solve(*fourierdata);
		pressure.solve(*fourierdata);
	}

	SingleParticle(const SingleParticle& particle) : data(particle.data), numNodes(particle.numNodes), numSeriesTerms(particle.numSeriesTerms) , flow(particle.numSeriesTerms), pressure(particle.numSeriesTerms), traction(&flow, &pressure)
	{
		data = new SphereData(*particle.data);
		fourierdata = new VSHSeries(*particle.fourierdata);
		*data = *particle.data;

		flow.solve(*fourierdata);
		pressure.solve(*fourierdata);
	}

	SingleParticle operator +(const SingleParticle & particle)
	{
		SingleParticle temp;

		for (int t = 0; t < NUMGLNODES; t++)
			for (int p = 0; p < NUMTRAPNODES; p++)
				temp.data[0][t][p] = data[0][t][p] + particle.data[0][t][p];

		*temp.fourierdata = *fourierdata + *particle.fourierdata;

		temp.numNodes = numNodes;

		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(*temp.fourierdata);
		temp.pressure.solve(*temp.fourierdata);
		temp.traction = StokesTractionHarmonic(&temp.flow, &temp.pressure);

		return temp;

	}

	SingleParticle operator -(const SingleParticle & particle)
	{
		SingleParticle temp;

		for (int t = 0; t < NUMGLNODES; t++)
			for (int p = 0; p < NUMTRAPNODES; p++)
				temp.data[0][t][p] = data[0][t][p] - particle.data[0][t][p];

		*temp.fourierdata = *fourierdata - *particle.fourierdata;

		temp.numNodes = numNodes;

		temp.numSeriesTerms = numSeriesTerms;

		temp.flow.solve(*temp.fourierdata);
		temp.pressure.solve(*temp.fourierdata);
		temp.traction = StokesTractionHarmonic(&temp.flow, &temp.pressure);

		return temp;

	}

	SingleParticle operator +=(const SingleParticle & particle)
	{
		*this = *this + particle;
		return *this;
	}

	SingleParticle operator *=(const double & a)
	{
		*this = (*this) * a;
		return *this;
	}

	SingleParticle operator *(const double& a)
	{
		SphereData* temp = new SphereData();
		for (int t = 0; t < NUMGLNODES; t++)
			for (int p = 0; p < NUMTRAPNODES; p++)
				temp[0][t][p] = a * data[0][t][p];

		return SingleParticle(temp,center , numSeriesTerms, numNodes  );
	}

	SingleParticle operator =(const SingleParticle & particle)
	{
		delete[] data;
		delete fourierdata;
		data = new SphereData(*particle.data);
		numNodes = particle.numNodes;
		numSeriesTerms = particle.numSeriesTerms;
		fourierdata =new VSHSeries( *particle.fourierdata);
		flow = SingleLayerHarmonic(*fourierdata);
		pressure = StokesPressure(*fourierdata);
		traction = StokesTractionHarmonic(&flow, &pressure);

		return *this;

	}

	double norm()
	{
		return sqrt(dot(*this));
	}

	double dot(SingleParticle& particle)
	{
		return L2InnerProductDiscrete(*data , *particle.data , numNodes);
	}

	static double dot(SingleParticle &particle, SingleParticle& qarticle)
	{
		return particle.dot(qarticle);
	}

	int getN()
	{
		return NUMGLNODES * NUMTRAPNODES * 3;
	}

	void axpy(SingleParticle* particle, double scale)
	{
		*this = *this + (*particle) * scale;
	}

	void operator ~()
	{
		delete[] data;
		delete[] fourierdata;
	}

	void refreshData()
	{
		delete fourierdata;
		fourierdata = new VSHSeries(numSeriesTerms, numNodes, center);
		fourierdata->approximate(*data , numSeriesTerms, numNodes);

		flow.solve(fourierdata[0]);
		pressure.solve(fourierdata[0]);
	}

};

RectCoord rotationCoefficient(const SphereData& data , RectCoord center, double radius)
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
	SingleParticle operator *(SingleParticle & particle)
	{
		SphereData* temp = new SphereData();

		RectCoord average = 0.25 / PI * IntegrateDiscrete(particle.data[0] , particle.numNodes);
		RectCoord rcoef = 0.25 / PI * rotationCoefficient(particle.data[0], particle.center, 1.0);

		for(int t = 0; t < NUMGLNODES; t++)
			for (int p = 0; p < particle.numNodes; p++)
			{
				SurfaceCoord s(PI / 2.0 * (GLnodes[t] + 1), 2.0 * PI * (double)p / (double)particle.numNodes);
				SphereCoord x(s);

				temp[0][t][p] = 0.5 * particle.data[0][t][p] + particle.traction(x) + average + cross(rcoef , SphereToRect(x) - x.center);

			}



		return SingleParticle(temp , particle.center, particle.numSeriesTerms , particle.numNodes);
	}


};

class RHSOperator
{
public:
	SingleParticle operator *( SingleParticle & particle)
	{
		SphereData* temp = new SphereData;


		for (int t = 0; t < NUMGLNODES; t++)
			for (int p = 0; p < NUMTRAPNODES; p++)
			{
				SurfaceCoord s(PI / 2.0 * (GLnodes[t] + 1), 2.0 * PI * (double)p / (double)NUMTRAPNODES);
				SphereCoord x(s);

				(*temp)[t][p] =-1.0 *( 0.5 * particle.data[0][t][p] + particle.traction(x));

			}


		return SingleParticle(temp,particle.center, particle.numSeriesTerms, particle.numNodes);
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

class Mtx
{
private:
	std::array<std::array<double, MAXITER>, MAXITER> data;
public:
	double& operator()(int m, int n)
	{
		return data[m][n];
	}
};

class ForceBalance : public SphericalVectorField
{
private:
	RectCoord force;
	RectCoord torque;
public:

	ForceBalance(RectCoord f , RectCoord t) : force(f){}

	RectCoord operator()(SphereCoord x)
	{
		return 0.25 / PI * (force + cross(torque , SphereToRect(x) - x.center));
	}
};



SingleParticle* SolveMobilitySingleSphere(RectCoord F , RectCoord T , int seriesSize)
{

	SphericalVectorField initialguess(&e_r);
		
	SingleParticle* solution = new SingleParticle(discretize(&initialguess) , RectCoord(), seriesSize , NUMTRAPNODES);
	LHSOperator* L = new LHSOperator();
	ForceBalance rho(F,T);
	SphereData* r;
	r = discretize(&rho);

	SingleParticle* rh = new SingleParticle(*r, RectCoord(), seriesSize, NUMTRAPNODES);
	RHSOperator R;
	SingleParticle* rhs = new SingleParticle();
	*rhs = R * (*rh);

	IdentityPreconditioner* I = new IdentityPreconditioner();

	GMRES(L, solution, rhs, I, 20, 5, 1e-10);


	
	SingleParticle* soln = new SingleParticle();
	*soln = *solution + *rh;

	delete L;
	delete I;
	delete solution;
	delete r;
    delete rh;
	delete rhs;
	return soln;

}