#ifdef __cplusplus
extern "C" {
#endif

struct plasma {
  double IonDensityPedistal;
  double IonDensitySeperatrix;
  double IonDensityOrigin;
  double IonTemperaturePedistal;
  double IonTemperatureSeperatrix;
  double IonTemperatureOrigin;
  double PedistalRadius;
  double IonDensityPeaking;
  double IonTemperaturePeaking;
  double MinorRadius;
  double MajorRadius;
  std::string PlasmaType;
  int    PlasmaId;
};

struct xs_params {
  double c[7];
};

void setup_(const double *ion_density_ped, const double *ion_density_sep,
	    const double *ion_density_origin, const double *ion_temp_ped,
	    const double *ion_temp_sep, const double *ion_temp_origin, 
	    const double *pedistal_rad, const double *ion_density_peak,
	    const double *ion_temp_peak, const double *minor_radius, 
	    const double *major_radius, const double *elongation, 
	    const double *triangularity, const double *shafranov, 
	    const char* plasma_type, const int *plasma_id,
	    const int* number_of_bins);

/*
 * Function to setup the plasma source in the first case.
 */
void setup_plasma_source(int number_of_bins,
			 const plasma plasma_parameters
 			 );
/*
 * function to calculate the ion density at a specific minor 
 * radius
 */
double ion_density(const plasma plasma_parameters,
		   double SampleRadius);

/*
 * function to calculate the ion temperature at a specific minor 
 * radius
 */
double ion_temperature(const plasma plasma_parameters,
		   double SampleRadius);

/*
 * function to determine the value of the dt xs cross sections at 
 * a specific ion temp
 */
double dt_xs(double ion_temp);

/*
 * sample the source, returns the minor radius sampled
 * expects new rn_store every call
 */
void sample_source_(const int *num_of_bins, const double *bin_width,double &sampled_radius, 
		    int &sampled_bin, double *rn_store1, double *rn_store2);

/*
 * sample the neutron energy  in MeV
 */
void sample_energy_(const int *bin_number, double *random_number1, double *random_number2,
		    double &energy_neutron);

/*
 * take the sampled minor radius and convert to cylindrical coordinates
 */ 
void convert_rad_to_rz_(const double *major_radius, const double *minor_radius,
			const double *elongation, const double *delta, 
			const double *shafranov, const double *minor_sampled,
			const double *rn_store, double &radius, double &height);

/*
 * convert partial cylindrical coords to xyz
 */
void convert_rz_to_xyz_(const double *r, const double *rn_store, 
		       double &x, double &y, double &toroidal_angle);

}
