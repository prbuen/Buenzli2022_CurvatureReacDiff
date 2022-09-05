/****************************************************************************** 
Main: model simulation

Author: Pascal R Buenzli, 2022

Description:
Buenzli PR and Simpson MJ (2022) Curvature dependences of wave propagation in reaction-diffusion models;
Preprint available at https://arxiv.org/abs/2112.00928

Requirements:
- scid library 0.3.x, https://code.dlang.org/packages/scid
- D compiler: ldc2 or dmd, see https://wiki.dlang.org/LDC
- Visualisation script written for Gnuplot 5.4.1
- Tested on MacOSX and Linux (Kubuntu 21.10)

Compilation:
From root directory:
	ldc2 src/*.d -od=build
or (optimised):
	ldc2 src/*.d -od=build -O3 -release -mcpu=native -ffast-math

Execution:
	./main
=> data created in `datadir/*.dat` (filenames are suffixed with time frame number)

Visualisation:
Edit plot.py (e.g. frame=70) and feed to gnuplot
	./plot.py | gnuplot
	pdflatex frame-0070
=> frame-0070.pdf
******************************************************************************/
import std.stdio;
import std.math;
import std.conv;
import std.range;
import std.exception;
import scid.matrix;

import utils: ModelInputsWrap, Data, Functions, curvature_centred_2, Vector;
import model_inputs;

alias MaterialProp = MatrixView!(Real); // Using scid's matrix library
alias Profile = Vector!(Real);

class State {
public:
	MaterialProp u; // 2D discretisation of cell density
	MatrixView!(bool) outside; //1 (true) outside the domain, 0 (false) inside
	static if (ModelInputs.calculate_curvature) MaterialProp kappa;
	static if (ModelInputs.calculate_normal_velocity) MaterialProp v;
private:
	MaterialProp u_old; // used for time stepping
}


// Model class
class ReacDiff2D {
protected:
	Real t=0.; // current time
	size_t iter=0; // current step() iteration
	size_t frame_iter=0; // current write_state() iteration
	State S; // current state
	ModelInputsWrap!(ModelInputs) p; // model parameters and input functions
	Data data; // data file directory and name patterns.

public:
	this()
	{
		p = new ModelInputsWrap!(ModelInputs)();		
		S = new State;
		
	   // read and assign parameter values
		p.readjust_discretisation_params(); // unnecessary?
		enforce(p.firstdatafile != "", std.conv.text("*** Error: no first data file provided."));
		stderr.writeln("#### Model: firstdatafile is '", p.firstdatafile, "'");

		// prepare output data file structure
		data = Data(p.firstdatafile);
		data.erase_datadir(); // start from a clean slate
		
		with(p)
		{			
			// output parameters to stderr and keep them as metadata in a datafile
			stderr.writeln(params() ~ "\n");
			stderr.writeln("#### Data written to = ", data.datadir, " ####");
			
			// p.metadata ~= params("# "); // comment out all
			p.metadata ~= params("# ", ["description", "model_version", "lambda"] ~ Functions); // comment out strings and functions, and parameter lambda (=python keyword)
			auto metadata = data.open("metadata.dat");
			metadata.write(p.metadata);
			metadata.flush();

			// ===== prepare state variables =====
			// allocate dynamic arrays
			S.u = matrix!Real(Nx, Ny);
			S.u_old = matrix!Real(Nx, Ny);
			S.outside = matrix!bool(Nx, Ny);

			static if(ModelInputs.calculate_curvature)
				S.kappa = matrix!Real(Nx, Ny);
			static if(ModelInputs.calculate_normal_velocity)
				S.v = matrix!Real(Nx, Ny);

			// initialise values
			foreach(i; 0 .. Nx)
				foreach(j; 0 .. Ny)
				{
					S.u[i,j] = p.u0(i,j,0);
					S.outside[i,j] = p.outside_domain(i,j,0);

					static if(ModelInputs.calculate_curvature)
					{
						Real Dx_u, Dy_u;
						S.kappa[i,j] = - curvature_centred_2(S.u, i, j, Dx_u, Dy_u)/dx;
					}
				}			

			// ===== write initial state =====
			write_state();			
		}		
	}

	// Repeatedly call `step()` until it returns false (stopping criterion)
	void run()
	{
		while(step()) {}
	}	
	
	/*
	   Write current `State` to datafiles
	*/
	void write_state()
	{
		auto u = data.open("u", frame_iter);
		assert(u.isOpen);

		// Potentially reduce spatial resolution; use binary float for fast data read/write
		foreach(i; iota(0, p.Nx, p.output_every_nth_dx))
			foreach(j; iota(0, p.Ny, p.output_every_nth_dy))
				u.writef("%r", to!float(S.u[i,j]));
		
		// domain mask
		if (frame_iter==0)
		{
			auto o = data.open("outside.dat");
			foreach(i; iota(0, p.Nx, p.output_every_nth_dx))
				foreach(j; iota(0, p.Ny, p.output_every_nth_dy))
					o.writef("%r", to!float(S.outside[i,j]));
		}
		
		static if (ModelInputs.calculate_curvature)
		{
			auto curv = data.open("curvature", frame_iter);
			foreach(i; iota(0, p.Nx, p.output_every_nth_dx))
				foreach(j; iota(0, p.Ny, p.output_every_nth_dy))
					curv.writef("%r", to!float(S.kappa[i,j]));
		}

		static if (ModelInputs.calculate_normal_velocity)
		{
			auto v = data.open("v", frame_iter);
			foreach(i; iota(0, p.Nx, p.output_every_nth_dx))
				foreach(j; iota(0, p.Ny, p.output_every_nth_dy))
					v.writef("%r", to!float(S.v[i,j]));
		}

		frame_iter++;
	}


	
	/*
	   Update the current `State` of the model
	*/
	bool step()
	{
		// keep previous values:
		with(p)
			with(S)
				foreach(i; 0 .. Nx)
					foreach(j; 0 .. Ny)
						u_old[i,j] = u[i,j];

		// update state:
		with(p)
			with (S)
			{
				foreach(i; 0 .. Nx)
					foreach(j; 0 .. Ny)
					{
						static if(ModelInputs.noflux_bc)
						{
							// no-flux boundary conditions (equivalent to Neumann BC here):
							// Applied using the ghost node method, i.e., for flux D(u)u_x, this gives
							// D(u_s) (u_{s+1}-u_{s-1}) / (2dx) = 0 => u_{s+1}=u_{s-1}, where u_{s} = boundary value.
							auto u_old_im1_j = (i==0?  u_old[i+1, j] : u_old[i-1,j]); // u_{-1,j} = u_{1,j} (i=0 boundary)
							auto u_old_ip1_j = (i==Nx-1? u_old[i-1,j] : u_old[i+1,j]); // u_{Nx,j} = u_{Nx-2,j} (i=Nx-1 boundary)
							
							auto u_old_i_jm1 = (j==0? u_old[i,j+1] : u_old[i, j-1]); // u_{i,-1} = u_{i,1} (j=0 boundary)
							auto u_old_i_jp1 = (j==Ny-1? u_old[i,j-1] : u_old[i,j+1]); // u_{i,Ny} = u_{i,Ny-2} (j=Ny-1 boundary)
						}
						else static if (ModelInputs.dirichlet_bc || ModelInputs.timedep_dirichlet_bc)
						{
							/* Dirichlet for arbitrary pores */
							// Set boundary values to u_initial
							outside[i,j] = outside_domain(i,j,t);
							
							static if (ModelInputs.noinnerBC)
								outside[i,j] = abs(u0(i,j,0)-u_initial) < 1e-10 && (i*dx-Lx/2.)^^2 + (j*dy-Ly/2.)^^2 > 1.1*R^^2;

							static if (ModelInputs.outward)
								outside[i,j] = outside[i,j] || i==0 || j==0 || i==Nx-1 || j==Ny-1;
								
							auto u_old_im1_j = (outside[i,j]? u0(i,j,t) : u_old[i-1,j]); // u(i-1,j) or boundary value
							auto u_old_ip1_j = (outside[i,j]? u0(i,j,t) : u_old[i+1,j]); // u(i+1,j) or boundary value
							
							auto u_old_i_jm1 = (outside[i,j]? u0(i,j,t) : u_old[i, j-1]); // u(i,j-1) or boundary value
							auto u_old_i_jp1 = (outside[i,j]? u0(i,j,t) : u_old[i,j+1]); // u(i,j+1) or boundary value
						}
						else
							assert(0);

						static if(ModelInputs.calculate_curvature || ModelInputs.calculate_normal_velocity)
						{
							Real Dx_u, Dy_u;
							kappa[i,j] = outside[i,j] || u[i,j] < p.u_min_for_curvature_velocity ? Real.nan : - curvature_centred_2(u, i, j, Dx_u, Dy_u)/dx;
							static if(ModelInputs.outward)
								kappa[i,j] = outside[i,j] || u[i,j] < p.u_min_for_curvature_velocity ? Real.nan : 1./sqrt((i*dx - Lx/2.)^^2 + (j*dy - Ly/2.)^^2);
						}

						uint diffstopfactor = 1;
						static if (ModelInputs.diffusion_stop)
							diffstopfactor = t >= t_diffusstop? 0 : 1;
						
						auto Dt_u = 1/(2.*dx^^2)*
							diffstopfactor*( (Diffus(u_old_ip1_j)+Diffus(u_old[i,j]))*(u_old_ip1_j-u_old[i,j]) -
							  (Diffus(u_old[i,j])+Diffus(u_old_im1_j))*(u_old[i,j]-u_old_im1_j) )
							
							+ 1/(2.*dy^^2)*
							diffstopfactor*( (Diffus(u_old_i_jp1)+Diffus(u_old[i,j]))*(u_old_i_jp1-u_old[i,j]) -
							  (Diffus(u_old[i,j])+Diffus(u_old_i_jm1))*(u_old[i,j]-u_old_i_jm1) )
							
							+ Reac(u_old[i,j]);

						static if(ModelInputs.calculate_normal_velocity)
							v[i,j] = outside[i,j] || u[i,j] < p.u_min_for_curvature_velocity? Real.nan : dx*Dt_u/sqrt(Dx_u^^2 + Dy_u^^2);
					
						u[i,j] += dt*Dt_u;
				}
								
				
			}

		
		// update time
		t += p.dt;
		iter++;
		if(iter % p.output_every_nth_dt == 0)
			write_state();
		
		return(t < p.t_end);		
	}
}


void main(string[] args)
{
	auto model = new ReacDiff2D;
	model.run();
}

