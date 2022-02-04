/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./custom.h"
#include <cmath>

void create_cell_types( void )
{
	// set the random seed
	SeedRandom( );
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = drag_update_cell_velocity;
	cell_defaults.phenotype.motility.lateral_restriction = parameters.doubles("lateral_restriction");
	cell_defaults.phenotype.motility.vertical_restriction = parameters.doubles("vertical_restriction");
	cell_defaults.phenotype.motility.forward_bias = parameters.doubles("forward_bias");

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 
	
	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
	cell_defaults.functions.update_phenotype = phenotype_function; 
	cell_defaults.functions.custom_cell_rule = custom_function; 
	cell_defaults.functions.contact_function = contact_function; 
	
	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 
	display_cell_definitions( std::cout ); 
	
	return; 
}

double locomotive_force_generator( )
{
	// random number generator to define cell velocities
	// based on anempirically obtained velocity distribution

	double random_value, force_value;
	double sigma = parameters.doubles("sigma");
	
	random_value = UniformRandom();
	force_value = sigma * pow(-2 * log(random_value), 0.5);

	return force_value;
}


void drag_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt )
{

	// find location of variables and base parameter values
	int ECM_density_index = microenvironment.find_density_index( "ECM" );

	// sample ECM
	double ECM_density = pCell->nearest_density_vector()[ECM_density_index];
	double dyn_viscosity;

	// get viscosity based on concentration
	if(ECM_density == 2.5)
	{
		dyn_viscosity = 7.96;
	}
	else if(ECM_density == 4.0)
	{
		dyn_viscosity = 18.42;
	}
	else if(ECM_density == 6.0)
	{
		dyn_viscosity = 39.15;
	}

	// update the speed value
	pCell->phenotype.motility.migration_speed = locomotive_force_generator();

	// update velocity
	standard_update_cell_velocity(pCell, phenotype, dt);

	// include the 1/vu (1/ECM density) term to consider friction
	pCell->velocity /= dyn_viscosity;

	return;

}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{

	
	Cell* pC;

	Cell_Definition* pCD = cell_definitions_by_index[0]; 

	double Xmin = -400; 
	double Ymin = microenvironment.mesh.bounding_box[1]; 
	double Zmin = -75; 

	double Xmax = 400; 
	double Ymax = microenvironment.mesh.bounding_box[4]; 
	double Zmax = 75; 
	

	double cell_radius = pCD->phenotype.geometry.radius; 
	double spacing = 2.0 * cell_radius * 2.0; 
	double half_space = 0.5*spacing; 
	double y_offset = sqrt(3.0)*half_space; 
	
	
	double x = Xmin + cell_radius; 
	double y = Ymin + cell_radius; 
	double z = Zmin + cell_radius; 
	
	int n = 0; 
	
	while( z <= Zmax - cell_radius )
	{
		while( x <= Xmax - cell_radius )
		{
			Cell* pC = create_cell( *pCD ); 
			pC->assign_position( x,y,z ); 
			
			x += spacing; 
		}
		x = Xmin + half_space; 
		
		n++; 
		z += y_offset; 
		if( n % 2 == 1 )
		{ x += half_space; }
	}
	
	
	// load cells from your CSV file (if enabled)
	// load_cells_from_pugixml(); 	
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; } 

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; } 