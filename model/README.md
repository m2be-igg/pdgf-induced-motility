# Model setup

The computational model was designed with PhysiCell (version 1.8.0), and the previously developed ECM extension.

In addition, the model was extended to provide more control over the migration patterns in 3D, by restricting movement in the vertical component and favouring forward/backwards movement.

## PhysiCell_phentype.cpp

New components were added to the motility class to enable **vertical** and **lateral restriction**. Cells can now also favour movement towards or against a given angle, defined as **forward_angle**. The probability of moving towards this angle is given by **forward_bias**.

```cpp
Motility::Motility()
{
	...

	vertical_restriction = 0.0;
	lateral_restriction = 0.0;
	forward_bias = 0.0;
	forward_angle = 3.1415926535897932384626433832795*0.5;

    ...
}
```

## PhysiCell_cell.cpp

The new components are used to **update the motility vector** at time intervals given by the persistence time.

```cpp
// Update the motility vector with the extended 3D movement function
void Cell::update_motility_vector( double dt_ )
{
    // Define a null vector if the cell is not moving
	if( phenotype.motility.is_motile == false )
	{
		phenotype.motility.motility_vector.assign( 3, 0.0 ); 
		return; 
	}
	

	if( UniformRandom() < dt_ / phenotype.motility.persistence_time || phenotype.motility.persistence_time < dt_ )
	{
		// If the update_bias_vector function is set, use it
		if( functions.update_migration_bias )
		{
			functions.update_migration_bias( this,phenotype,dt_ );
			phenotype.motility.forward_angle = atan2(phenotype.motility.migration_bias_direction[1], phenotype.motility.migration_bias_direction[0]);
		}

		// Limit cell motility in the x component
		double temp_angle;

		if( UniformRandom()  < phenotype.motility.forward_bias ) {
			temp_angle = 3.1415926535897932384626433832795*((1 - phenotype.motility.lateral_restriction)*UniformRandom()
					+ (-0.5 + phenotype.motility.lateral_restriction/2)) + phenotype.motility.forward_angle;
		} else {
			temp_angle = 3.1415926535897932384626433832795*((-1 + phenotype.motility.lateral_restriction)*UniformRandom()
					+ (1.5 - phenotype.motility.lateral_restriction/2)) + phenotype.motility.forward_angle;
		}
		
		// Limit cell motility in the z component
		double temp_phi = 3.1415926535897932384626433832795*((1 - phenotype.motility.vertical_restriction)*UniformRandom()
				+ phenotype.motility.vertical_restriction/2);

		double sin_phi = sin(temp_phi);
		double cos_phi = cos(temp_phi);
		
		if( phenotype.motility.restrict_to_2D == true )
		{ 
			sin_phi = 1.0; 
			cos_phi = 0.0;
		}
		
		std::vector<double> randvec; 
		randvec.resize(3,sin_phi); 
		
		randvec[0] *= cos( temp_angle ); // cos(theta)*sin(phi)
		randvec[1] *= sin( temp_angle ); // sin(theta)*sin(phi)
		randvec[2] = cos_phi; //  cos(phi)

		phenotype.motility.motility_vector = randvec;

		normalize( &(phenotype.motility.motility_vector) );

		phenotype.motility.motility_vector *= phenotype.motility.migration_speed;
	}	
	return; 
} 
```