/*
Model inputs: parameters and given functions

Author: Pascal R Buenzli, 2022

Description:
Several parameter sets are defined. These correspond to different application of the governing equations.
*/
module model_inputs;

import std.math;
import std.algorithm;
import std.conv;
import std.exception;
import scid.matrix;

import utils: DefSwitch, adjust_to_dx_coarse_given_dx, adjust_to_dx_coarse_given_length;

alias Real = double;


// ===== Uncomment the model to compile =====
/* Square 1x1 pore */
// version = Model_0_N; // no-flux BC
version = Model_0_D; // Dirichlet BC

/* Circular pore, Dirichlet */
// version = Model_2; // inward
// version = Model_2_out;
// version = Model_2_out_noinnerBC; // outward, no inner boundary, only initial condition
// version = Model_2_Skellam_D; // inward, time-dependent Dirichlet BC
// version = Model_2_AllenCahn; // Allen-Cahn with u\in [-1,1] instead of [0,1]
// version = Model_2_Allee; // strong Allee with u \in [0,1]

/* Square pore of side a, Dirichlet */
// version = Model_3; // see Model_0
// version = Model_3_out;

/* Star-like pore (n branches, e.g. n=3), Dirichlet */
// version = Model_4;
// version = Model_4_out;

/* Cross-like pore, Dirichlet*/
// version = Model_6;
// version = Model_6_out;

// ===== Global options =====
// version = noflux_bc;
// version = dirichlet_bc;
// version = timedep_dirichlet_bc;
version = calculate_curvature;
version = calculate_normal_velocity;
// version = diffusion_stop; // stop diffusion after some time
// version = reac_zero; // no source/sink term (diffusion equation)
// version = reac_linear; // linear source term F(u)=u (Skellam)
version = reac_logistic; // logistic source term F(u) = u(1-u)
// version = reac_u2_1mu; // F(u) = u^2(1-u)
// version = reac_1mu; // F(u) = 1-u (i.e. linearised around F(1)=0, t->inf) (saturated growth)
// version = reac_u_uma_1mu; // F(u) = u(u-a)(1-u) (strong Allee effect)
// version = outward; // expanding out instead of infilling
// version = noinnerBC; // for outward fronts, don't implement any BC in the centre, only IC



/********************************************************************************************/
// automatically set model parameters according to model version chosen above:
version(Model_0_N) { alias ModelInputs = ModelInputs_0; version=noflux_bc; };
version(Model_0_D) { alias ModelInputs = ModelInputs_0; version=dirichlet_bc; };
version(Model_2) { alias ModelInputs = ModelInputs_2; version=dirichlet_bc; };
version(Model_2_Skellam_D) { alias ModelInputs = ModelInputs_2_Skellam; version=timedep_dirichlet_bc; version=reac_1mu;};
version(Model_2_AllenCahn) { alias ModelInputs = ModelInputs_2_AllenCahn; version=dirichlet_bc; version=reac_u_1mu2;};
version(Model_2_Allee) { alias ModelInputs = ModelInputs_2_Allee; version=dirichlet_bc; version=reac_u_uma_1mu;};
version(Model_3) { alias ModelInputs = ModelInputs_3; version=dirichlet_bc; };
version(Model_4) { alias ModelInputs = ModelInputs_4; version=dirichlet_bc; };
version(Model_6) { alias ModelInputs = ModelInputs_6; version=dirichlet_bc; };
version(Model_2_out) { alias ModelInputs = ModelInputs_2_out; version=dirichlet_bc; version=outward;};
version(Model_2_out_noinnerBC) { alias ModelInputs = ModelInputs_2_out; version=dirichlet_bc; version=outward; version=noinnerBC;};
version(Model_3_out) { alias ModelInputs = ModelInputs_3_out; version=dirichlet_bc; version=outward;};
version(Model_4_out) { alias ModelInputs = ModelInputs_4_out; version=dirichlet_bc; version=outward;};
version(Model_6_out) { alias ModelInputs = ModelInputs_6_out; version=dirichlet_bc; version=outward;};
/********************************************************************************************/
	

class ModelInputs_0 {
public:
	// version parameters
	string metadata; // will be filled by params
	string firstdatafile;

	// global options
	// export compilation switches for use in other modules (DefSwitch is from pb.config)
	mixin(DefSwitch!("noflux_bc"));
	mixin(DefSwitch!("dirichlet_bc"));
	mixin(DefSwitch!("timedep_dirichlet_bc"));
	mixin(DefSwitch!("diffusion_stop")); // stop diffusion after t_diffusstop

	mixin(DefSwitch!("reac_zero"));
	mixin(DefSwitch!("reac_linear"));
	mixin(DefSwitch!("reac_logistic"));
	mixin(DefSwitch!("reac_u2_1mu"));
	mixin(DefSwitch!("reac_1mu"));
	mixin(DefSwitch!("reac_u_1mu2"));
	mixin(DefSwitch!("reac_u_uma_1mu"));
	
	mixin(DefSwitch!("calculate_curvature"));
	mixin(DefSwitch!("calculate_normal_velocity"));

	mixin(DefSwitch!("outward"));
	mixin(DefSwitch!("noinnerBC"));

	// time discretisation:
	Real t_end;
	Real dt; // computation time step
	double dt_frame; // visualisation time step

	// space discretisation:
	// We assume periodic boundary condition (PBC) along x and will store Nx values: x_0, ..., x_{N-1}, where x_N==x_0	
	Real Lx; // domain length
	Real Ly; // domain height	
	Real dx;
	Real dy;
	Real dx_coarse;
	Real dy_coarse;

	// model parameters
	Real u_initial; // initial density at the edge; u* in the paper
	Real u_c; // threshold density for red contour
	Real D; // diffusivity
	Real R; // ModelInputs_2: radius; ModelInputs_4: baseline radius
	Real A; // Allee threshold in ModelInputs_2_Allee, F(u)=u(u-A)(1-u), A<1
	Real a; // ModelInputs_3: edge length of initial square pore, must be < 1; ModelInputs_4: size of perturbation, must have R+a < Lx, Ly; ModelInputs_6: scaling factor
	uint n; // ModelInputs_4: number of branches in the 'star'
	Real t_diffusstop; // time at which diffusion stops
	Real u_min_for_curvature_velocity; // density threshold below which no curvature nor normal velocity is computed (low res)
	
	// dependent simulation parameters
    size_t output_every_nth_dt, Nt_frame, Nt;
	size_t output_every_nth_dx, output_every_nth_dy, Nx, Ny, Nx_coarse, Ny_coarse;
	
	/* To be called each time dx, dt etc. are changed */
	void readjust_discretisation_params()
	{
		// calculates output_every_nth_dt, Nt_coarse, Nt, realigns dt_coarse and t_end:
		adjust_to_dx_coarse_given_dx(t_end, dt, Nt, dt_frame, Nt_frame, output_every_nth_dt);

		// calculates output_every_nth_dx, Nx_coarse, Nx, realigns dx_coarse and x_max with PBC.
		adjust_to_dx_coarse_given_length(Lx, dx, Nx, dx_coarse, Nx_coarse, output_every_nth_dx);
			
		// calculates output_every_nth_dy, Ny_coarse, Ny, realigns dx_coarse and y_max with PBC.
		adjust_to_dx_coarse_given_length(Ly, dy, Ny, dy_coarse, Ny_coarse, output_every_nth_dy);
	}	

	this()
	{
		// sanity checks
		static if (noflux_bc && dirichlet_bc) enforce (0, "noflux_bc and dirichlet_bc can't be set jointly");
		static if (! (noflux_bc || dirichlet_bc || timedep_dirichlet_bc)) enforce (0, "noflux_bc or (time-dep or const) dirichlet_bc must be set");
		static if (int(reac_zero) + int(reac_linear) + int(reac_logistic) + int(reac_u2_1mu) + int(reac_1mu) + int(reac_u_1mu2) + int(reac_u_uma_1mu) != 1)
			enforce(0, std.conv.text("Can't have no or several reaction terms: reac_zero=", reac_zero, "; reac_linear=", reac_linear, "; reac_logistic=", reac_logistic, "; reac_u2_1mu=", reac_u2_1mu, "; reac_1mu=", reac_1mu, "; reac_u_1mu2=", reac_u_1mu2, "; reac_u_uma_1mu=", reac_u_uma_1mu));
		
		/*****  workaround ******/
		// non-NaN values are needed here even if not used
		R = 1e10; // Only used >= ModelInputs_2
		A = 1e10; // Only used in ModelInputs_2_Allee
		a = 1e10; // >= ModelInputs_3; Ensure there is a layer of boundary sites with u_initial
		n = 0; // ModelInputs_4;
		t_diffusstop = 9.0;
		/***********************/

		// parameter values
		firstdatafile = "datadir/u_0.dat";
		
		D = 0.005;
		t_end = 20.;
		dt = 0.001;
		dt_frame = 0.1;

		Lx = 2.;
		Ly = 2.;
		dx = 0.01;
		dy = 0.01;
		dx_coarse = 0.01;
		dy_coarse = 0.01;
		u_min_for_curvature_velocity = 0.0;

		u_c = 0.5;
		u_initial = 1;
		R = 1.0;
	}

	Real u0(size_t i, size_t j, Real t)
	{
		static if (noflux_bc)
		{
			if (i<=1 || i>=Nx-2 || j<=1 || j>=Ny-2)
				return u_initial;
			else
				return 0;
		}
		else static if (dirichlet_bc)
		{
			if (i==0 || i==Nx-1 || j==0 || j==Ny-1)
				return u_initial;
			else
				return 0;
		}
		else
			assert(0);
	}

	bool outside_domain(size_t i, size_t j, Real t)
	{
		return i==0 || i==Nx-1 || j==0 || j==Ny-1;
	}
	
	// Diffusivity
	Real Diffus(Real u)
	{
		return D;
	}

	Real Reac(Real u)
	{

		static if (reac_zero)
			return 0;
		else static if (reac_linear)
			return u;
		else static if (reac_logistic)
			return u*(1.0-u);
		else static if (reac_u2_1mu)
			return u^^2*(1.0-u);
		else static if (reac_1mu)
			return 1.0-u;
		else static if (reac_u_1mu2)
			return u*(1.0-u^^2);
		else static if (reac_u_uma_1mu)
			return u*(u-A)*(1.0-u)/(1.0-A); // div by 1-A so when u->1, F(u)~1-u
		else
			enforce(0, "No other reaction terms implemented.");
		return 0;
	}
}


class ModelInputs_2 : ModelInputs_0 {
public:
	this() { }

	override Real u0(size_t i, size_t j, Real t)
	{
		static if (noflux_bc)
		{
			enforce(0, "Model_2 (circular pore) with Neumann BC: not implemented in 2D");
			if (i<=1 || i>=Nx-2 || j<=1 || j>=Ny-2)
				return u_initial;
			else
				return 0;
		}
		else static if (dirichlet_bc)
		{
			if ( (i*dx - Lx/2.)^^2 + (j*dy - Ly/2.)^^2 >= R^^2 )
				return u_initial;
			else
				return 0;			
		}
		else static if (timedep_dirichlet_bc)
		{
			if ( (i*dx - Lx/2.)^^2 + (j*dy - Ly/2.)^^2 >= R^^2 )
				return u_initial*exp(t);
			else
				return 0;			
		}
		else
		{
			enforce(0);
			return 0;
		}
	}
}

class ModelInputs_2_Skellam : ModelInputs_2 {
public:
	this()
	{
		D = 0.05;
		t_end = 9.0;
		dt = 0.0001;
		dt_frame = 0.03;
	}
}


class ModelInputs_2_out : ModelInputs_2 {
public:
	this()
	{
		D = 0.02;
		t_end = 20.;
		dt = 0.002;
		dt_frame = 0.1;

		Lx = 16.;
		Ly = 16.;
		dx = 0.05;
		dy = 0.05;
		dx_coarse = 0.05;
		dy_coarse = 0.05;
		u_min_for_curvature_velocity = 1e-6;
	}

	override Real u0(size_t i, size_t j, Real t)
	{
		static if (noflux_bc)
		{
			enforce(0, "Model_2 (circular pore) with Neumann BC: not implemented");
			if (i<=1 || i>=Nx-2 || j<=1 || j>=Ny-2)
				return u_initial;
			else
				return 0;
		}
		else static if (dirichlet_bc)
		{
			if ( (i*dx - Lx/2.)^^2 + (j*dy - Ly/2.)^^2 <= R^^2 )
				return u_initial;
			else
				return 0;			
		}
		else
		{
			enforce(0);
			return 0;
		}
	}

}


class ModelInputs_2_AllenCahn : ModelInputs_0 {
public:
	this()
	{
	}

	override Real u0(size_t i, size_t j, Real t)
	{
		// wave-pinning:
		static if (dirichlet_bc)
		{
			auto r2 = (i*dx - Lx/2.)^^2 + (j*dy - Ly/2.)^^2;
			if ( r2 >= R^^2 )
				return u_initial; // Dirichlet bc
			else
				return 1-2*cos(PI*sqrt(r2)/(2*R));
		}		
		else
		{
			enforce(0, "Model_2_AllenCahn (circular pore) with non-Dirichlet BC not implemented");
			return 0;
		}
	}

	override bool outside_domain(size_t i, size_t j, Real t)
	{
		if ( (i*dx - Lx/2.)^^2 + (j*dy - Ly/2.)^^2 >= R^^2 )
			return 1;
		else
			return 0;
	}
}


class ModelInputs_2_Allee : ModelInputs_0 {
public:
	this()
	{		
		t_end = 60.;
		A = 2./5.; // Allee threshold
	}

	override Real u0(size_t i, size_t j, Real t)
	{
		static if (dirichlet_bc)
		{
			if ( (i*dx - Lx/2.)^^2 + (j*dy - Ly/2.)^^2 >= R^^2 )
				return u_initial; // Dirichlet bc
			else
				return 0;
		}
		else
		{
			enforce(0, "Model_2_Allee (circular pore) with other BC: not implemented");
			return 0;
		}
	}

	override bool outside_domain(size_t i, size_t j, Real t)
	{
		if ( (i*dx - Lx/2.)^^2 + (j*dy - Ly/2.)^^2 >= R^^2 )
			return 1;
		else
			return 0;
	}
}


class ModelInputs_3 : ModelInputs_0 {
public:
	this()
	{		
		t_end = 15.;
		u_min_for_curvature_velocity = 1e-6;
		a=1.98; //edge length of initial square pore, must be < Lx, e.g. Lx-2*dx to ensure there is a layer of boundary sites with u_initial. Otherwise see Model0
	}

	override Real u0(size_t i, size_t j, Real t)
	{
		static if (noflux_bc)
		{
			enforce(0, "Model_3 (inner square pore) with Neumann BC: not implemented");
		}
		else static if (dirichlet_bc)
		{
			if ( abs(i*dx - Lx/2.) >= a/2. || abs(j*dy - Ly/2.) >= a/2. )
				return u_initial;
			else
				return 0;			
		}
		else
			enforce(0);
		return 0;
	}

	override bool outside_domain(size_t i, size_t j, Real t)
	{
			if ( abs(i*dx - Lx/2.) >= a/2. || abs(j*dy - Ly/2.) >= a/2. )
			return 1;
		else
			return 0;
	}	
}

class ModelInputs_3_out : ModelInputs_3 {
public:
	// Real a; // edge length of initial square pore, must be < 1

	this()
	{
		D = 0.02;
		t_end = 20.;
		dt = 0.002;
		dt_frame = 0.1;

		Lx = 16.;
		Ly = 16.;
		dx = 0.05;
		dy = 0.05;
		dx_coarse = 0.05;
		dy_coarse = 0.05;
		u_min_for_curvature_velocity = 1e-6;
	}

	override Real u0(size_t i, size_t j, Real t)
	{
		static if (noflux_bc)
		{
			enforce(0, "Model_3 (inner square pore) with Neumann BC: not implemented");
		}
		else static if (dirichlet_bc)
		{
			if ( abs(i*dx - Lx/2.) <= a/2. && abs(j*dy - Ly/2.) <= a/2. )
				return u_initial;
			else
				return 0;			
		}
		else
			enforce(0);
		return 0;
	}

	override bool outside_domain(size_t i, size_t j, Real t)
	{
		return super.outside_domain(i,j,t)==true? false : true;
	}
}

class ModelInputs_4 : ModelInputs_0 {
public:
	this()
	{
		t_end=15.;
		R = 0.65; // baseline radius
		a = 0.35; // size of perturbation. We must have R + a < Lx,Ly
		n = 3; // number of branches
	}

	override Real u0(size_t i, size_t j, Real t)
	{
		static if (noflux_bc)
		{
			enforce(0, "Model_4 (3 petal) with Neumann BC not implemented");
		}
		else static if (dirichlet_bc)
		{
			auto x = i*dx - Lx/2; // x position relative to pore centre
			auto y = j*dx - Ly/2; // y position relative to pore centre
			auto norm2 = x^^2 + y^^2;
			if( norm2 <= (R - a)^^2)
				return 0; // avoid finding the polar angle theta at the origin
			else
			{
				auto theta = atan2(y,x);  // polar angle of (x,y)
				if (norm2 < (R + a*cos(n*theta))^^2)
					return 0;
				else
					return u_initial;
			}
		}
		else
			enforce(0);
		return 0;	
	}

	override bool outside_domain(size_t i, size_t j, Real t)
	{
			auto x = i*dx - Lx/2; // x position relative to pore centre
			auto y = j*dx - Ly/2; // y position relative to pore centre
			auto norm2 = x^^2 + y^^2;
			if( norm2 <= (R - a)^^2)
				return 0; // avoid finding the polar angle theta at the origin
			else
			{
				auto theta = atan2(y,x);  // polar angle of (x,y)
				if (norm2 < (R + a*cos(n*theta))^^2)
					return 0;
				else
					return 1;
			}
	}
}


class ModelInputs_4_out : ModelInputs_4 {
public:
	// Real R; // baseline radius
	// Real a; // size of perturbation. We must have R+a < Lx and Ly
	// uint n; // number of branches

	this()
	{
		D = 0.02;
		t_end = 20.;
		dt = 0.002;
		dt_frame = 0.1;
		Lx = 16.;
		Ly = 16.;
		dx = 0.05;
		dy = 0.05;
		dx_coarse = 0.05;
		dy_coarse = 0.05;
		u_min_for_curvature_velocity = 1e-6;
	}

	override Real u0(size_t i, size_t j, Real t)
	{
		static if (noflux_bc)
		{
			enforce(0, "Model_4_out (3 petal) with Neumann BC: not implemented");
		}
		else static if (dirichlet_bc)
		{
			auto x = i*dx - Lx/2; // x position relative to pore centre
			auto y = j*dx - Ly/2; // y position relative to pore centre
			auto norm2 = x^^2 + y^^2;
			if( norm2 <= (R - a)^^2)
				return u_initial; // avoid finding the polar angle theta at the origin
			else
			{
				auto theta = atan2(y,x);  // polar angle of (x,y)
				if (norm2 > (R + a*cos(n*theta))^^2)
					return 0;
				else
					return u_initial;
			}
		}
		else
			enforce(0);
		return 0;	
	}


	override bool outside_domain(size_t i, size_t j, Real t)
	{
		return super.outside_domain(i,j,t)==true? false : true;
	}
}



class ModelInputs_6 : ModelInputs_0 {
public:
	this()
	{
		t_end = 7.;
		a = 1.0; // scaling factor
	}

	override Real u0(size_t i, size_t j, Real t)
	{
		static if (noflux_bc)
		{
			enforce(0, "Model_6 (inner square pore) with Neumann BC: not implemented");
		}
		else static if (dirichlet_bc)
		{
			auto x = i*dx - Lx/2; // x position relative to pore centre
			auto y = j*dx - Ly/2; // y position relative to pore centre
			auto norm2 = x^^2 + y^^2;
			if( norm2 <= sqrt(2.)*dx)
				return 0; // avoid finding the polar angle theta at the origin
			else
			{
				auto theta = atan2(y,x);  // polar angle of (x,y)
				// The cross curve is R(theta)=sin^4(theta) + cos^4(theta). It fits in the square [-1,1]^2. Above we have already offset the origin to the pore centre.
				auto R_theta = sin(theta)^^4 + cos(theta)^^4;
				
				if (norm2 < R_theta^^2)
					return 0;
				else
					return u_initial;
			}
		}
		else
	
			enforce(0);
	
		return 0;
	}


	override bool outside_domain(size_t i, size_t j, Real t)
	{
		auto x = i*dx - Lx/2; // x position relative to pore centre
			auto y = j*dx - Ly/2; // y position relative to pore centre
			auto norm2 = x^^2 + y^^2;
			if( norm2 <= sqrt(2.)*dx)
				return 0; // avoid finding the polar angle theta at the origin
			else
			{
				auto theta = atan2(y,x);  // polar angle of (x,y)
				// The cross curve is R(theta)=sin^4(theta) + cos^4(theta). It fits in the square [-1,1]^2. Above we have already offset the origin to the pore centre.
				auto R_theta = sin(theta)^^4 + cos(theta)^^4;
				
				if (norm2 < R_theta^^2)
					return 0;
				else
					return 1;
			}
	}
}


class ModelInputs_6_out : ModelInputs_6 {
public:
	// Real a; // scaling factor

	this()
	{
		D = 0.02;
		t_end = 20.;
		dt = 0.002;
		dt_frame = 0.1;
		Lx = 16.;
		Ly = 16.;
		dx = 0.05;
		dy = 0.05;
		dx_coarse = 0.05;
		dy_coarse = 0.05;
		u_min_for_curvature_velocity = 1e-6;
	}

	override Real u0(size_t i, size_t j, Real t)
	{
		static if (noflux_bc)
		{
			enforce(0, "Model_6 (inner square pore) with Neumann BC: not implemented");
		}
		else static if (dirichlet_bc)
		{
			auto x = i*dx - Lx/2; // x position relative to pore centre
			auto y = j*dx - Ly/2; // y position relative to pore centre
			auto norm2 = x^^2 + y^^2;
			if( norm2 <= sqrt(2.)*dx)
				return u_initial; // avoid problems finding the polar angle theta at the origin
			else
			{
				auto theta = atan2(y,x);  // polar angle of (x,y), between -pi and pi, well behaved at x=0
				// The cross curve is R(theta)=sin^4(theta) + cos^4(theta). It fits in the square [-1,1]^2. Above we have already offset the origin to the pore centre.
				auto R_theta = sin(theta)^^4 + cos(theta)^^4;
				
				if (norm2 > R_theta^^2)
					return 0;
				else
					return u_initial;
			}
		}
		else
	
			enforce(0);
	
		return 0;
	}

	override bool outside_domain(size_t i, size_t j, Real t)
	{
		return super.outside_domain(i,j,t)==true? false : true;
	}	
}
