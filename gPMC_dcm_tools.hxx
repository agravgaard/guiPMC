#include "gdcmReader.h"
#include "gdcmAttribute.h"

#ifdef GPMC_UI_H
#include <QString>
#include <QDir>
#include "cl.hpp"
#include "range_modulator_data.hxx"
#include "itkEuler3DTransform.h"
#include <chrono>
#include <random>
#include "gPMC_commandline.hxx"
// DCM writer stuff
#include "itkGDCMImageIO.h"
#include "itkGDCMSeriesFileNames.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageSeriesReader.h"
#include "itkImageSeriesWriter.h"
#include "itkResampleImageFilter.h"
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
#include <itksys/SystemTools.hxx>
#include "gdcmUIDGenerator.h"
#include <string>
#include <sstream>
#endif

#define M_PI 3.1415926535897932384626433832795028841971693993751
#define N 1000000

//                        N points      x array       y array      x point      y point
int point_in_aperture(int nvert, float *vertx, float *verty, float testx, float testy)
{
	int i, j, c = 0;
	for (i = 0, j = nvert - 1; i < nvert; j = i++) {
		if (((verty[i] > testy) != (verty[j] > testy)) &&
			(testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]))
			c = !c;
	}
	return c; // 0=outside, 1=inside
}

//           range modulator angles     range modulator heights    material     lead or alu angles         lead or alu heights
//std::tuple<std::vector<const double>, std::vector<const double>, const char*, std::vector<const double>, std::vector<const double>>
// args: (j, n_lead_alu_only, mod_start, std::get<3>(range_modulator)[j], std::get<3>(range_modulator)[j + 1], std::get<0>(range_modulator)[0])
bool is_lead_alu_angle_in_modulated_angle_space(const size_t cur_idx, const size_t n_lead_alu_only, const double mod_start,
	const double cur_lead_alu_angle, const double nxt_lead_alu_angle, const double first_rng_mod_angle)
{
	if (cur_idx < (n_lead_alu_only - 1)){
		if (mod_start <= cur_lead_alu_angle){
			if (first_rng_mod_angle >= nxt_lead_alu_angle)
				return true;
			else if (first_rng_mod_angle >= cur_lead_alu_angle)
				return true;
		}
		else if ((first_rng_mod_angle >= nxt_lead_alu_angle) && (mod_start < nxt_lead_alu_angle))
			return true;
	}
	else if ((first_rng_mod_angle >= cur_lead_alu_angle) && (mod_start <= cur_lead_alu_angle))
		return true;

	return false;
}
// args: (j, n_rng_mod, mod_start, mod_end, std::get<0>(range_modulator)[j_RM], std::get<0>(range_modulator)[j_RM + 1])
bool is_carbon_angle_in_modulated_angle_space(const size_t cur_idx, const size_t n_rng_mod, const double mod_start, const double mod_end,
	const double cur_rng_mod_angle, const double nxt_rng_mod_angle)
{
	if (cur_idx < (n_rng_mod - 1)){
		if (mod_start <= cur_rng_mod_angle){
			if (mod_end >= nxt_rng_mod_angle)
				return true;
			else if (mod_end >= cur_rng_mod_angle)
				return true;
		}
		else if ((mod_end >= nxt_rng_mod_angle) && (mod_start < nxt_rng_mod_angle))
			return true;
	}
	else if ((mod_start <= cur_rng_mod_angle) && (mod_end >= cur_rng_mod_angle))
		return true;

	return false;
}

template <size_t NperSpot>
size_t GetNFromDicom(const std::string dicom_path){
	// check dicom integrity
	if (dicom_path.empty()) return 0;

	gdcm::Reader reader;
	reader.SetFileName(dicom_path.c_str());
	if (!reader.Read())
	{
		std::cout << "Reading dicom plan failed!" << std::endl;
		return 69;
	}
	gdcm::File &file = reader.GetFile();
	gdcm::DataSet &ds = file.GetDataSet();

	size_t total_n_spots = 0;
	const gdcm::DataElement &beam_seq_tag = ds.GetDataElement(gdcm::Tag(0x300a, 0x3a2));

	gdcm::SmartPointer<gdcm::SequenceOfItems> beam_seq = beam_seq_tag.GetValueAsSQ();
	for (unsigned int i = 0; i < beam_seq->GetNumberOfItems(); ++i){
		const gdcm::Item &it_beams = beam_seq->GetItem(i + 1);
		gdcm::Attribute< 0x300a, 0x308> at_scan_mode;
		at_scan_mode.SetFromDataElement(it_beams.GetDataElement(at_scan_mode.GetTag()));
		const char* scan_mode = at_scan_mode.GetValue();

		if (!strcmp(scan_mode, "NONE")){ // passive scatter
			// Range compensators:
			const gdcm::DataElement &rng_comp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x2ea));
			gdcm::SmartPointer<gdcm::SequenceOfItems> rng_comp_seq = rng_comp_seq_tag.GetValueAsSQ();
			const gdcm::Item &it_rng_comp = rng_comp_seq->GetItem(1);

			gdcm::Attribute< 0x300a, 0xe7> at_n_rows; // number of rows in rng comp. -> x-direction.
			at_n_rows.SetFromDataElement(it_rng_comp.GetDataElement(at_n_rows.GetTag()));
			const unsigned int n_rows = (unsigned)at_n_rows.GetValue(); // value can't be signed, ignore warning
			gdcm::Attribute< 0x300a, 0xe8> at_n_cols; // number of columns in rng comp. -> y-direction.
			at_n_cols.SetFromDataElement(it_rng_comp.GetDataElement(at_n_cols.GetTag()));
			const unsigned int n_cols = (unsigned)at_n_cols.GetValue(); // value can't be signed, ignore warning

			gdcm::Attribute< 0x300a, 0xe9> at_spacing; // spacing of rows and columns in rng comp.
			at_spacing.SetFromDataElement(it_rng_comp.GetDataElement(at_spacing.GetTag()));
			const double* spacing = at_spacing.GetValues(); // mm
			gdcm::Attribute< 0x300a, 0xea> at_offset; // position of 0,0 in rng comp. -> upper lefthand.
			at_offset.SetFromDataElement(it_rng_comp.GetDataElement(at_offset.GetTag()));
			const double* offset = at_offset.GetValues(); // mm

			// Block:
			const gdcm::DataElement &block_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a6));
			gdcm::SmartPointer<gdcm::SequenceOfItems> block_seq = block_seq_tag.GetValueAsSQ();
			const gdcm::Item &it_block = block_seq->GetItem(1);

			gdcm::Attribute< 0x300a, 0x104> at_block_n_points;
			at_block_n_points.SetFromDataElement(it_block.GetDataElement(at_block_n_points.GetTag()));
			const int block_n_points = at_block_n_points.GetValue(); // number of x-y pairs defining block
			assert(block_n_points > 0);

			gdcm::Attribute< 0x300a, 0x106> at_block_points;
			at_block_points.SetFromDataElement(it_block.GetDataElement(at_block_points.GetTag()));
			const double* block_points = at_block_points.GetValues(); // x-y pairs defining block
			float* block_points_x = new float[(unsigned)block_n_points + 1];
			float* block_points_y = new float[(unsigned)block_n_points + 1];

#pragma omp parallel for
			for (int j = 0; j < block_n_points; j++){
				block_points_x[j] = block_points[j * 2];
				block_points_y[j] = block_points[j * 2 + 1];
			}
			block_points_x[(unsigned)block_n_points] = block_points[0]; //close the loop
			block_points_y[(unsigned)block_n_points] = block_points[1];

			size_t n_spots_in_aperture = 0;
			for (size_t j = 0; j < n_cols; j++){
				const cl_float cur_col_pos = j * spacing[1] - offset[1];
				for (size_t k = 0; k < n_rows; k++){ // k + j*n_rows
					if (1 == point_in_aperture(block_n_points, block_points_x, block_points_y,
						k * spacing[0] + offset[0], cur_col_pos))
						n_spots_in_aperture++;
				}
			}

			delete[] block_points_x;
			delete[] block_points_y;
			// end block

			size_t n_mod_wepl = 0;
			const gdcm::DataElement &cp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a8));
			gdcm::SmartPointer<gdcm::SequenceOfItems> control_point_seq = cp_seq_tag.GetValueAsSQ();

			const gdcm::Item &it_control_points = control_point_seq->GetItem(1);

			// Range modulator settings:
			const gdcm::DataElement &rng_mod_set_seq_tag = it_control_points.GetDataElement(gdcm::Tag(0x300a, 0x380));
			gdcm::SmartPointer<gdcm::SequenceOfItems> rng_mod_set_seq = rng_mod_set_seq_tag.GetValueAsSQ();

			for (unsigned int k = 0; k < rng_mod_set_seq->GetNumberOfItems(); ++k){
				const gdcm::Item &it_rng_mod_set = rng_mod_set_seq->GetItem(k + 1);
				gdcm::Attribute< 0x300a, 0x382> at_mod_start; // range modulator start index
				at_mod_start.SetFromDataElement(it_rng_mod_set.GetDataElement(at_mod_start.GetTag()));
				const unsigned int mod_start = at_mod_start.GetValue();

				gdcm::Attribute< 0x300a, 0x384> at_mod_end; // range modulator end index
				at_mod_end.SetFromDataElement(it_rng_mod_set.GetDataElement(at_mod_end.GetTag()));
				const unsigned int mod_end = at_mod_end.GetValue();

				// Range modulator:
				const gdcm::DataElement &rng_mod_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x342));
				gdcm::SmartPointer<gdcm::SequenceOfItems> rng_mod_seq = rng_mod_seq_tag.GetValueAsSQ();

				const gdcm::Item &it_rng_mod = rng_mod_seq->GetItem(1);

				gdcm::Attribute< 0x300a, 0x346> at_mod_ID; // range modulator type
				at_mod_ID.SetFromDataElement(it_rng_mod.GetDataElement(at_mod_ID.GetTag()));
				const char* mod_ID = at_mod_ID.GetValue();
				if (!isdigit(mod_ID[3])){
					std::cout << "\a" << "Non-standard format of range modulator id: " << mod_ID << " source code edit may be necessary." << std::endl;
					return 0;
				}

				const size_t ID_num = (size_t)(mod_ID[3] - '0'); // apparently the best way to convert from char to int

				//         range modulator angles    range modulator heights   material     lead or alu angles        lead or alu heights
				std::tuple<std::vector<const double>, std::vector<const double>, const char*, std::vector<const double>, std::vector<const double>>
					range_modulator = get_range_modulator(ID_num);

				size_t n_lead_alu_only = 0;
				for (size_t j = 0; j < std::get<3>(range_modulator).size(); j++)
					if (std::get<3>(range_modulator)[j] < std::get<0>(range_modulator)[0] && std::get<3>(range_modulator)[j] > mod_start)
						n_lead_alu_only++;

				const size_t n_rng_mod = std::get<0>(range_modulator).size() + n_lead_alu_only;

				size_t n_non_zero_weights = 0;
				// Range modulator energy difference
				const double first_rng_mod_angle = std::get<0>(range_modulator)[0];
				for (size_t j = 0; j < n_lead_alu_only; j++){
					const double nxt_angle = (j < (n_lead_alu_only - 1)) ? std::get<3>(range_modulator)[j + 1] : 0.0;
					if (is_lead_alu_angle_in_modulated_angle_space(j, n_lead_alu_only, mod_start, std::get<3>(range_modulator)[j],
						nxt_angle, first_rng_mod_angle))
						n_non_zero_weights++;
				}

				for (size_t j = n_lead_alu_only; j < n_rng_mod; j++){
					const size_t j_RM = j - n_lead_alu_only;
					const double nxt_angle = (j < (n_rng_mod - 1)) ? std::get<0>(range_modulator)[j_RM + 1] : 0.0;
					if (is_carbon_angle_in_modulated_angle_space(j, n_rng_mod, mod_start, mod_end, std::get<0>(range_modulator)[j_RM], nxt_angle))
						n_non_zero_weights++;
				}

				n_mod_wepl += (n_non_zero_weights * NperSpot);
			}
			total_n_spots += n_spots_in_aperture * n_mod_wepl;
		}
		else // spot scanning
		{
			const gdcm::DataElement &cp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a8));
			gdcm::SmartPointer<gdcm::SequenceOfItems> control_point_seq = cp_seq_tag.GetValueAsSQ();

			for (unsigned int j = 0; j < control_point_seq->GetNumberOfItems(); ++j){
				const gdcm::Item &it_control_points = control_point_seq->GetItem(j + 1);
				gdcm::Attribute< 0x300a, 0x392> at_n_spots;
				at_n_spots.SetFromDataElement(it_control_points.GetDataElement(at_n_spots.GetTag()));
				total_n_spots += at_n_spots.GetValue()  * NperSpot;
			}
		}
	}
	return total_n_spots;
}

double inverse_cdf(double random){
	for (size_t i = 0; i < (sizeof(cdf_xs) / sizeof(double)); i++)
		if (random <= cdf_xs[i]) // cdf_xs goes from 0 to 1
			return cdf_val[i];
	return cdf_val[sizeof(cdf_val) / sizeof(double) - 1]; // if we have not returned yet, return the last element (impossible in theory, but compiler doesn't know that)
}

template <size_t NperSpot>
size_t initFromSpotScanning(gdcm::SmartPointer<gdcm::SequenceOfItems> beam_seq,
	cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight, size_t& total_spots){
	typedef itk::Euler3DTransform< double > TransformType;
	TransformType::ParametersType fixedParam(3); //rotation center
	fixedParam.put(0, 0);
	fixedParam.put(1, 0);
	fixedParam.put(2, 0);

	const double halfC = M_PI / 180.0;

	for (unsigned int i = 0; i < beam_seq->GetNumberOfItems(); ++i) { // LOOPING BEAMS, assumes there are at least one beam
		const gdcm::Item &it_beams = beam_seq->GetItem(i + 1);
		gdcm::Attribute<0x300a, 0xC6> at_rt_type;
		at_rt_type.SetFromDataElement(it_beams.GetDataElement(at_rt_type.GetTag()));
		const std::string rt_type(at_rt_type.GetValue());

		float ion_modifier = 1.0; // Energy is float
		if (rt_type.find("ION") != std::string::npos)
		{
			gdcm::Attribute<0x300a, 0x306> at_Z; //charge(306) or atom#(304) better?
			at_Z.SetFromDataElement(it_beams.GetDataElement(at_Z.GetTag()));
			gdcm::Attribute<0x300a, 0x302> at_A;
			at_A.SetFromDataElement(it_beams.GetDataElement(at_A.GetTag()));
			ion_modifier = float(at_Z.GetValue()) / float(at_A.GetValue());
			std::cout << "Warning: not proton plan! -> multiplying by Z/A for approximation: " << ion_modifier << std::endl;
		}

		const gdcm::DataElement &cp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a8));
		gdcm::SmartPointer<gdcm::SequenceOfItems> control_point_seq = cp_seq_tag.GetValueAsSQ();

		const gdcm::Item &first_control_points = control_point_seq->GetItem(1);
		// While we're on the first control point we get the data that may only be defined here
		gdcm::Attribute<0x300a, 0x11e> at_gantry;
		at_gantry.SetFromDataElement(first_control_points.GetDataElement(at_gantry.GetTag()));
		const double gantry = at_gantry.GetValue() * halfC;

		gdcm::Attribute<0x300a, 0x122> at_couch;
		at_couch.SetFromDataElement(first_control_points.GetDataElement(at_couch.GetTag()));
		const double couch = at_couch.GetValue() * halfC;

		gdcm::Attribute<0x300a, 0x12c> at_isocenter;
		at_isocenter.SetFromDataElement(first_control_points.GetDataElement(at_isocenter.GetTag()));
		const __m128 isocenter = {
			at_isocenter.GetValues()[0],
			-at_isocenter.GetValues()[1], // notice this minus! it migth be different for other plans than the syngo ones.
			at_isocenter.GetValues()[2],
			0.0f
		};

		// RAnge Shifter
		const gdcm::DataElement &range_shifter_seq_tag = first_control_points.GetDataElement(gdcm::Tag(0x300a, 0x360));
		gdcm::SmartPointer<gdcm::SequenceOfItems> range_shifter_seq = range_shifter_seq_tag.GetValueAsSQ();
		const gdcm::Item &it_rng_shifter = range_shifter_seq->GetItem(1);
		gdcm::Attribute<0x300a, 0x364> at_sid; // Isocenter to Rangeshifter Distance
		at_sid.SetFromDataElement(it_rng_shifter.GetDataElement(at_sid.GetTag()));
		const __m128 sid3 = _mm_set1_ps(at_sid.GetValue()); // mm //cl_float3 to use vector operations
		// const float sid = at_sid.GetValue();

		gdcm::Attribute<0x300a, 0x366> at_rswet; // Rangeshifter Water Equivalent Thickness
		at_rswet.SetFromDataElement(it_rng_shifter.GetDataElement(at_rswet.GetTag()));
		const float rs_wet = at_rswet.GetValue();

		const __m128 direction = {
			std::sin(gantry)*std::cos(couch),
			std::cos(gantry),
			std::sin(gantry)*std::sin(couch),
			0.0f
		};

		//            ( x )                              ( x' )       ( x )             [a b c]
		// IF point = ( y ) THEN point after rotation is ( y' ) = A * ( y ) , WHERE A = [d e f]
		//            ( z )                              ( z' )       ( z )             [g h i]
		TransformType::Pointer transform = TransformType::New();
		//transform->SetRotation(0, couch, gantry);
		//      ASSUMING:         3., 2., 1.
		transform->SetRotation(couch, 0, gantry);
		transform->SetFixedParameters(fixedParam); //Center of the Transform
		transform->SetComputeZYX(true); // Just to be sure it behaves as expected
		const itk::Matrix<double, 3U, 3U> A = transform->GetMatrix();

		const __m128 beam_offset = _mm_sub_ps(isocenter, _mm_mul_ps(direction, sid3));
		/*const cl_float3 beam_offset = {
			isocenter[0] - direction.x * sid,
			-isocenter[1] - direction.y * sid,
			isocenter[2] - direction.z * sid
			};*/

		for (unsigned int j = 0; j < control_point_seq->GetNumberOfItems(); ++j) { // LOOPING-CONTROL POINTS, assumes there are at least one control point per beam.
			const gdcm::Item &it_control_points = control_point_seq->GetItem(j + 1);
			gdcm::Attribute<0x300a, 0x114> at_nom_beam_energy;
			at_nom_beam_energy.SetFromDataElement(it_control_points.GetDataElement(at_nom_beam_energy.GetTag()));
			const float nom_beam_energy = at_nom_beam_energy.GetValue() * ion_modifier - rs_wet * bethe(at_nom_beam_energy.GetValue() * ion_modifier); // gPMC uses float and it's not used for anything else.

			gdcm::Attribute< 0x300a, 0x392> at_n_spots;
			at_n_spots.SetFromDataElement(it_control_points.GetDataElement(at_n_spots.GetTag()));
			const int n_spots = at_n_spots.GetValue();
			assert(n_spots >= 0);

			gdcm::Attribute<0x300a, 0x394> at_spot_pos_map; // x, y - positions of spots
			at_spot_pos_map.SetFromDataElement(it_control_points.GetDataElement(at_spot_pos_map.GetTag()));
			const float* p_spot_pos_map = at_spot_pos_map.GetValues(); // gPMC uses float and it's not used for anything else.
			const std::vector<float> spot_pos_map{ p_spot_pos_map, p_spot_pos_map + at_spot_pos_map.GetNumberOfValues() };

			gdcm::Attribute<0x300a, 0x396> at_spot_meterset_w;
			at_spot_meterset_w.SetFromDataElement(it_control_points.GetDataElement(at_spot_meterset_w.GetTag()));
			const float* p_spot_w = at_spot_meterset_w.GetValues(); // gPMC uses float and it's not used for anything else.
			const std::vector<float> spot_meterset_w{ p_spot_w, p_spot_w + at_spot_meterset_w.GetNumberOfValues() };

			gdcm::Attribute<0x300a, 0x398> at_spot_size;
			at_spot_size.SetFromDataElement(it_control_points.GetDataElement(at_spot_size.GetTag()));
			const float* spot_size = at_spot_size.GetValues();

#pragma omp parallel for // (Ignore signed/unsigned warning! omp works only with int not unsigned)
			for (int i_spot_xy = 0; i_spot_xy < (2 * n_spots); i_spot_xy += 2)
			{
				std::minstd_rand0 g1((unsigned)std::clock() + i_spot_xy);  // minstd_rand0 is a standard linear_congruential_engine
				const float g1_max = g1.max();
				const float rand_to_angle = M_PI / g1_max; // only half circle because inverse_cdf gives a value from a normal distribution with µ=0
				const float rand_to_radius_x = 0.1f * spot_size[0];
				const float rand_to_radius_y = 0.1f * spot_size[1];

				for (size_t k = 0; k < NperSpot; k++)
				{
					const unsigned int i_local = total_spots + (i_spot_xy / 2 * NperSpot) + k;
					// Dir(ection) is the opposite (unit)vector of that pointing from iso->source
					// Theta is the angle from z->phi and phi is angle from x->y
					// theta is gantry if z is towards anterior and phi is couch if y is cranial and x is lateral right
					float* tmp_float3 = new float[4];
					_mm_store_ps(tmp_float3, direction); // x-lateral   = -sin(theta)*cos(phi)
					// dir[i_local].s[1] = direction[1]; // y-AP        = -cos(theta)
					// dir[i_local].s[2] = direction[2]; // z-CC        = -sin(theta)*sin(phi)
					dir[i_local] = { tmp_float3[0], tmp_float3[1], tmp_float3[2] };

					const float rand_sin = sin(float(g1()) * rand_to_angle);
					const float rand_cos = sqrt(1 - rand_sin * rand_sin);

					const float rand_radius_x = inverse_cdf(float(g1()) / g1_max) * rand_to_radius_x;
					const float rand_radius_y = inverse_cdf(float(g1()) / g1_max) * rand_to_radius_y;

					const cl_float2 point = {
						spot_pos_map[(size_t)i_spot_xy] + rand_cos * rand_radius_x,
						spot_pos_map[(size_t)i_spot_xy + 1] + rand_sin * rand_radius_y
					};

					const __m128 trans_point = { //   y=0, so only z of point
						A(0, 0) * point.x + A(0, 2) * point.y,
						A(1, 0) * point.x + A(1, 2) * point.y,
						A(2, 0) * point.x + A(2, 2) * point.y,
						0.0f // just to be 4 floats wide
					};

					_mm_store_ps(tmp_float3, _mm_add_ps(trans_point, beam_offset));
					pos[i_local] = { tmp_float3[0], tmp_float3[1], tmp_float3[2] };

					T[i_local] = nom_beam_energy;
					weight[i_local] = spot_meterset_w[(size_t)i_spot_xy / 2];
				}
			} // end for and omp parallel for
			total_spots += n_spots * NperSpot; // for enabeling full parallelism
		}
	}

	return total_spots;
}

size_t spots_in_block(cl_float2* spot_map, const size_t n_rows, const size_t n_cols,
	const double* block_points, const int block_n_points,
	const double* spacing, const double* offset)
{
	assert(block_n_points > 0 && spacing > 0);
	float* block_points_x = new float[(unsigned)block_n_points + 1];
	float* block_points_y = new float[(unsigned)block_n_points + 1];

#pragma omp parallel for
	for (int j = 0; j < block_n_points; j++){
		block_points_x[j] = block_points[j * 2];
		block_points_y[j] = block_points[j * 2 + 1];
	}
	block_points_x[(unsigned)block_n_points] = block_points[0]; //close the loop
	block_points_y[(unsigned)block_n_points] = block_points[1];

	size_t n_spots_in_aperture = 0;
	for (size_t j = 0; j < n_cols; j++){
		const cl_float cur_col_pos = j * spacing[1] - offset[1];
		for (size_t k = 0; k < n_rows; k++){ // k + j*n_rows
			spot_map[n_spots_in_aperture].x = k * spacing[0] + offset[0]; //row pos // unsure if + or - offset
			spot_map[n_spots_in_aperture].y = cur_col_pos; //col pos
			if (1 == point_in_aperture(block_n_points, block_points_x, block_points_y,
				spot_map[n_spots_in_aperture].x, spot_map[n_spots_in_aperture].y))
				n_spots_in_aperture++;
		}
	}

	delete[] block_points_x;
	delete[] block_points_y;

	return n_spots_in_aperture;
}

std::tuple<size_t, double*, double*> range_modulator_wheel_to_stopping_power(
	const std::tuple<std::vector<const double>, std::vector<const double>, const char*, std::vector<const double>, std::vector<const double>> range_modulator,
	const size_t ID_num, const float nom_beam_energy, const double mod_start, const double mod_end)
{
	size_t j_W = 0;

	size_t n_lead_alu_only = 0;
	for (size_t j = 0; j < std::get<3>(range_modulator).size(); j++)
		if (std::get<3>(range_modulator)[j] < std::get<0>(range_modulator)[0] && std::get<3>(range_modulator)[j] > mod_start)
			n_lead_alu_only++;

	const size_t n_rng_mod = std::get<0>(range_modulator).size() + n_lead_alu_only;
	double* mod_energy_diff = new double[n_rng_mod];
	double* mod_energy_weigth = new double[n_rng_mod];

	// Range modulator energy difference
	for (size_t j = 0; j < n_lead_alu_only; j++){
		if (std::get<3>(range_modulator)[j] > mod_start){
			mod_energy_diff[j_W] = ((ID_num != 9) ? bethe_lead(nom_beam_energy) : bethe_aluminium(nom_beam_energy))
				* std::get<4>(range_modulator)[j];
		}
		if (j < (n_lead_alu_only - 1)){
			if (mod_start <= std::get<3>(range_modulator)[j]){
				if (std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j + 1])
					mod_energy_weigth[j_W] = (std::get<3>(range_modulator)[j + 1] - std::get<3>(range_modulator)[j]);
				else if (std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j])
					mod_energy_weigth[j_W] = (std::get<0>(range_modulator)[0] - std::get<3>(range_modulator)[j]);
				else
					mod_energy_weigth[j_W] = 0.0;
			}
			else if ((std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j + 1]) && (mod_start < std::get<3>(range_modulator)[j + 1]))
				mod_energy_weigth[j_W] = (std::get<3>(range_modulator)[j + 1] - mod_start);
			else
				mod_energy_weigth[j_W] = 0.0;
		}
		else if ((mod_start <= std::get<3>(range_modulator)[j]) && (std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j]))
			mod_energy_weigth[j_W] = (std::get<0>(range_modulator)[0] - std::get<3>(range_modulator)[j]);
		else
			mod_energy_weigth[j_W] = 0.0;
		if (mod_energy_weigth[j_W] != 0.0)
			j_W++;
	}

	for (size_t j = n_lead_alu_only; j < n_rng_mod; j++){
		size_t j_RM = j - n_lead_alu_only;
		if (std::get<1>(range_modulator)[j_RM] == NULL)
			mod_energy_diff[j_W] = -nom_beam_energy * 0.99; // kill this point with 1% uncertainty
		else if (!strcmp(std::get<2>(range_modulator), "carbon"))
			mod_energy_diff[j_W] = -bethe_carbon(nom_beam_energy) * std::get<1>(range_modulator)[j_RM];
		else if (!strcmp(std::get<2>(range_modulator), "carbon2"))
			mod_energy_diff[j_W] = -bethe_carbon2(nom_beam_energy) * std::get<1>(range_modulator)[j_RM];
		else if (!strcmp(std::get<2>(range_modulator), "lexan"))
			mod_energy_diff[j_W] = -bethe_lexan(nom_beam_energy) * std::get<1>(range_modulator)[j_RM];

		for (size_t k = 0; k < std::get<3>(range_modulator).size(); k++){ // angle lead or alu foil
			if (std::get<3>(range_modulator)[k] == std::get<0>(range_modulator)[j_RM]){
				mod_energy_diff[j_W] -= ((ID_num != 9) ? bethe_lead(nom_beam_energy) : bethe_aluminium(nom_beam_energy))
					* std::get<4>(range_modulator)[k];
			}
		}
		if (j < (n_rng_mod - 1)){
			if (mod_start <= std::get<0>(range_modulator)[j_RM]){
				if (mod_end >= std::get<0>(range_modulator)[j_RM + 1])
					mod_energy_weigth[j_W] = (std::get<0>(range_modulator)[j_RM + 1] - std::get<0>(range_modulator)[j_RM]);
				else if (mod_end >= std::get<0>(range_modulator)[j_RM])
					mod_energy_weigth[j_W] = (mod_end - std::get<0>(range_modulator)[j_RM]);
				else
					mod_energy_weigth[j_W] = 0.0;
			}
			else if ((mod_end >= std::get<0>(range_modulator)[j_RM + 1]) && (mod_start < std::get<0>(range_modulator)[j_RM + 1]))
				mod_energy_weigth[j_W] = (std::get<0>(range_modulator)[j_RM + 1] - mod_start);
			else
				mod_energy_weigth[j_W] = 0.0;
		}
		else
		{
			if (mod_end < std::get<0>(range_modulator)[j_RM])
				mod_energy_weigth[j_W] = 0.0;
			else if (mod_start <= std::get<0>(range_modulator)[j_RM])
				mod_energy_weigth[j_W] = (mod_end - std::get<0>(range_modulator)[j_RM]);
		}
		if (mod_energy_weigth[j_W] != 0.0)
			j_W++;
	}
	return std::make_tuple(j_W, mod_energy_diff, mod_energy_weigth);
}

template <size_t NperSpot>
size_t initFromPassiveScatter(gdcm::SmartPointer<gdcm::SequenceOfItems> beam_seq,
	cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight, size_t& n_total_spots){
	typedef itk::Euler3DTransform< double > TransformType;
	TransformType::ParametersType fixedParam(3); //rotation center
	fixedParam.put(0, 0);
	fixedParam.put(1, 0);
	fixedParam.put(2, 0);

	const double halfC = M_PI / 180.0;
	for (unsigned int i = 0; i < beam_seq->GetNumberOfItems(); ++i) { // LOOPING BEAMS, assumes there are at least one beam
		const gdcm::Item &it_beams = beam_seq->GetItem(i + 1);
		gdcm::Attribute<0x300a, 0xC6> at_rt_type;
		at_rt_type.SetFromDataElement(it_beams.GetDataElement(at_rt_type.GetTag()));
		const std::string rt_type(at_rt_type.GetValue());

		float ion_modifier = 1.0; // Energy is float
		if (rt_type.find("ION") != std::string::npos)
		{
			gdcm::Attribute<0x300a, 0x306> at_Z; //charge(306) or atom#(304) better?
			at_Z.SetFromDataElement(it_beams.GetDataElement(at_Z.GetTag()));
			gdcm::Attribute<0x300a, 0x302> at_A;
			at_A.SetFromDataElement(it_beams.GetDataElement(at_A.GetTag()));
			ion_modifier = float(at_Z.GetValue()) / float(at_A.GetValue());
			std::cout << "Warning: not proton plan! -> multiplying by Z/A for approximation: " << ion_modifier << std::endl;
		}
		// Range compensators:
		const gdcm::DataElement &rng_comp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x2ea));
		gdcm::SmartPointer<gdcm::SequenceOfItems> rng_comp_seq = rng_comp_seq_tag.GetValueAsSQ();

		const gdcm::Item &it_rng_comp = rng_comp_seq->GetItem(1);

		gdcm::Attribute< 0x300a, 0x2e4> at_cid; // compensator isocenter distance in mm
		at_cid.SetFromDataElement(it_rng_comp.GetDataElement(at_cid.GetTag()));
		const __m128 sid3 = _mm_set1_ps(at_cid.GetValue()); // mm //cl_float3 to use vector operations
		// for source geometry

		gdcm::Attribute< 0x300a, 0x2e0> at_div; // compensator divergence (ABSENT or PRESENT)
		at_div.SetFromDataElement(it_rng_comp.GetDataElement(at_div.GetTag()));
		const char* div = at_div.GetValue();
		if (!strcmp(div, "PRESENT"))
			std::cout << "WARNING: COMPENSATOR DIVERGENCE PRESENT FOR BUT NOT ACCOUNTED FOR!" << std::endl;

		gdcm::Attribute< 0x300a, 0xe7> at_n_rows; // number of rows in rng comp. -> x-direction.
		at_n_rows.SetFromDataElement(it_rng_comp.GetDataElement(at_n_rows.GetTag()));
		const unsigned int n_rows = (unsigned)at_n_rows.GetValue(); // value can't be signed, ignore warning
		gdcm::Attribute< 0x300a, 0xe8> at_n_cols; // number of columns in rng comp. -> y-direction.
		at_n_cols.SetFromDataElement(it_rng_comp.GetDataElement(at_n_cols.GetTag()));
		const unsigned int n_cols = (unsigned)at_n_cols.GetValue(); // value can't be signed, ignore warning

		size_t n_spots_max = n_rows*n_cols;
		cl_float2* spot_map = new cl_float2[n_spots_max];

		gdcm::Attribute< 0x300a, 0xe9> at_spacing; // spacing of rows and columns in rng comp.
		at_spacing.SetFromDataElement(it_rng_comp.GetDataElement(at_spacing.GetTag()));
		const double* spacing = at_spacing.GetValues(); // mm
		gdcm::Attribute< 0x300a, 0xea> at_offset; // position of 0,0 in rng comp. -> upper lefthand.
		at_offset.SetFromDataElement(it_rng_comp.GetDataElement(at_offset.GetTag()));
		const double* offset = at_offset.GetValues(); // mm

		// Block:
		const gdcm::DataElement &block_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a6));
		gdcm::SmartPointer<gdcm::SequenceOfItems> block_seq = block_seq_tag.GetValueAsSQ();
		const gdcm::Item &it_block = block_seq->GetItem(1);

		gdcm::Attribute< 0x300a, 0xf8> at_block_type;
		at_block_type.SetFromDataElement(it_block.GetDataElement(at_block_type.GetTag()));
		const char* block_type = at_block_type.GetValue(); // APERTURE or SHIELDING
		if (!strcmp(block_type, "SHIELDING"))
			std::cout << "WARNING: BLOCK TAGGED " << block_type << ", ONLY APERTURE MODE IMPLEMENTED!" << std::endl;

		gdcm::Attribute< 0x300a, 0xe1> at_block_mat;
		at_block_mat.SetFromDataElement(it_block.GetDataElement(at_block_mat.GetTag()));
		const char* block_mat = at_block_mat.GetValue(); // Brass probably
		gdcm::Attribute< 0x300a, 0x100> at_block_thicc;
		at_block_thicc.SetFromDataElement(it_block.GetDataElement(at_block_thicc.GetTag()));
		const double block_thicc = at_block_thicc.GetValue(); // thickness in mm
		std::cout << "Block material: " << block_mat << " of " << block_thicc << " mm thickness";
		std::cout << " assumed to stop beam completely." << std::endl;

		gdcm::Attribute< 0x300a, 0x104> at_block_n_points;
		at_block_n_points.SetFromDataElement(it_block.GetDataElement(at_block_n_points.GetTag()));
		const int block_n_points = at_block_n_points.GetValue(); // number of x-y pairs defining block

		gdcm::Attribute< 0x300a, 0x106> at_block_points;
		at_block_points.SetFromDataElement(it_block.GetDataElement(at_block_points.GetTag()));
		const double* block_points = at_block_points.GetValues(); // x-y pairs defining block

		size_t n_spots_in_aperture = spots_in_block(spot_map, n_rows, n_cols, block_points, block_n_points, spacing, offset);
		//                                          ^passed by reference to generate spot_map

		gdcm::Attribute< 0x300a, 0x2e5> at_hex_offset; // column offset of in rng comp. -> only applicaple for hexagonal compensators.
		at_hex_offset.SetFromDataElement(it_rng_comp.GetDataElement(at_hex_offset.GetTag()));
		const double hex_offset = at_hex_offset.GetValue();
		if (hex_offset != 0)
			std::cout << "WARNING: HEXAGONAL COMPENSATOR DETECTED FOR BUT NOT IMPLEMENTED PROPERLY! expect errors." << std::endl;

		gdcm::Attribute< 0x300a, 0xec> at_thicc; // thickness of rng comp. in mm
		at_thicc.SetFromDataElement(it_rng_comp.GetDataElement(at_thicc.GetTag()));
		const double* thicc = at_thicc.GetValues();
		// DOUBLE SIDED triggers "isocenter to compensator distances" 300a,0x2e6

		gdcm::Attribute< 0x300a, 0x2e7> at_dedx; // Compensator Linear Stopping Power Ratio, ...
		// relative to water, at the beam energy specified by the Nominal Beam Energy (300A,0114) ...
		// of the first Control Point of the Ion Control Point Sequence (300A,03A8).
		at_dedx.SetFromDataElement(it_rng_comp.GetDataElement(at_dedx.GetTag()));
		const double dedx_rel = at_dedx.GetValue(); //relative to water

		gdcm::Attribute< 0x300a, 0xec> at_milling_diameter; // Compensator Milling Tool Diameter
		at_milling_diameter.SetFromDataElement(it_rng_comp.GetDataElement(at_milling_diameter.GetTag()));
		const double mill = at_milling_diameter.GetValue();

		// Range modulator:
		const gdcm::DataElement &rng_mod_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x342));
		gdcm::SmartPointer<gdcm::SequenceOfItems> rng_mod_seq = rng_mod_seq_tag.GetValueAsSQ();

		const gdcm::Item &it_rng_mod = rng_mod_seq->GetItem(1);

		gdcm::Attribute< 0x300a, 0x346> at_mod_ID; // range modulator type
		at_mod_ID.SetFromDataElement(it_rng_mod.GetDataElement(at_mod_ID.GetTag()));
		const char* mod_ID = at_mod_ID.GetValue();
		if (!isdigit(mod_ID[3])){
			std::cout << "\a" << "Non-standard format of range modulator id: " << mod_ID << " source code edit may be necessary." << std::endl;
			return 0;
		}

		const size_t ID_num = (size_t)(mod_ID[3] - '0'); // apparently the best way to convert from char to int (Ignore signed/unsigned warning!)

		//         range modulator angles    range modulator heights   material     lead or alu angles        lead or alu heights
		std::tuple<std::vector<const double>, std::vector<const double>, const char*, std::vector<const double>, std::vector<const double>>
			range_modulator = get_range_modulator(ID_num);

		gdcm::Attribute< 0x300a, 0x348> at_mod_type; // range modulator type
		at_mod_type.SetFromDataElement(it_rng_mod.GetDataElement(at_mod_type.GetTag()));
		const char* mod_type = at_mod_type.GetValue();
		if (strcmp(mod_type, "WHL_MODWEIGHTS"))
			std::cout << "\a" << "WARNING: WHL_MODWEIGHTS NOT DEFINED EXPECT FATAL ERRORS!" << std::endl;

		//control points:
		const gdcm::DataElement &cp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a8));
		gdcm::SmartPointer<gdcm::SequenceOfItems> control_point_seq = cp_seq_tag.GetValueAsSQ();
		const gdcm::Item &it_control_points = control_point_seq->GetItem(1);
		// While we're on the first control point we get the data that may only be defined here

		// Range modulator settings:
		const gdcm::DataElement &rng_mod_set_seq_tag = it_control_points.GetDataElement(gdcm::Tag(0x300a, 0x380));
		gdcm::SmartPointer<gdcm::SequenceOfItems> rng_mod_set_seq = rng_mod_set_seq_tag.GetValueAsSQ();

		const gdcm::Item &it_rng_mod_set = rng_mod_set_seq->GetItem(1);

		gdcm::Attribute< 0x300a, 0x382> at_mod_start; // range modulator start index
		at_mod_start.SetFromDataElement(it_rng_mod_set.GetDataElement(at_mod_start.GetTag()));
		const unsigned int mod_start = at_mod_start.GetValue();

		gdcm::Attribute< 0x300a, 0x384> at_mod_end; // range modulator end index
		at_mod_end.SetFromDataElement(it_rng_mod_set.GetDataElement(at_mod_end.GetTag()));
		const unsigned int mod_end = at_mod_end.GetValue();

		gdcm::Attribute<0x300a, 0x11e> at_gantry;
		at_gantry.SetFromDataElement(it_control_points.GetDataElement(at_gantry.GetTag()));
		const double gantry = at_gantry.GetValue() * halfC;

		gdcm::Attribute<0x300a, 0x122> at_couch;
		at_couch.SetFromDataElement(it_control_points.GetDataElement(at_couch.GetTag()));
		const double couch = at_couch.GetValue() * halfC;

		gdcm::Attribute<0x300a, 0x12c> at_isocenter;
		at_isocenter.SetFromDataElement(it_control_points.GetDataElement(at_isocenter.GetTag()));
		const __m128 isocenter = {
			at_isocenter.GetValues()[0],
			at_isocenter.GetValues()[1],
			at_isocenter.GetValues()[2],
			0.0f // just to be 4 floats wide
		};

		const __m128 direction = {
			std::sin(gantry)*std::cos(couch),
			std::cos(gantry),
			std::sin(gantry)*std::sin(couch),
			0.0f // just to be 4 floats wide
		};

		//            ( x )                              ( x' )       ( x )             [a b c]
		// IF point = ( y ) THEN point after rotation is ( y' ) = A * ( y ) , WHERE A = [d e f]
		//            ( z )                              ( z' )       ( z )             [g h i]
		TransformType::Pointer transform = TransformType::New();
		//transform->SetRotation(0, couch, gantry);
		//      ASSUMING:         3., 2., 1.
		transform->SetRotation(couch, 0, gantry);
		transform->SetFixedParameters(fixedParam); //Center of the Transform
		transform->SetComputeZYX(true); // Just to be sure it behaves as expected
		const itk::Matrix<double, 3U, 3U> A = transform->GetMatrix();

		// Beam offset [mm]
		const __m128 beam_offset = _mm_sub_ps(isocenter, _mm_mul_ps(direction, sid3));

		gdcm::Attribute<0x300a, 0x114> at_nom_beam_energy;
		at_nom_beam_energy.SetFromDataElement(it_control_points.GetDataElement(at_nom_beam_energy.GetTag()));
		const float nom_beam_energy = at_nom_beam_energy.GetValue() * ion_modifier; // gPMC uses float and it's not used for anything else.

		// Using geometry defined in range_modulator_data.hxx calculate the beam modulation of the beam. (assumes contant rotation of wheel)
		std::tuple<size_t, double*, double*> mod_wheel_tuple = range_modulator_wheel_to_stopping_power(
			range_modulator, ID_num, nom_beam_energy, mod_start, mod_end);

		size_t j_W = std::get<0>(mod_wheel_tuple);
		double* mod_energy_diff = std::get<1>(mod_wheel_tuple);
		double* mod_energy_weigth = std::get<2>(mod_wheel_tuple);

		// increment to second control point to get weight of beam:
		const gdcm::Item &second_control_point = control_point_seq->GetItem(2);

		gdcm::Attribute<0x300a, 0x134> at_weight;
		at_weight.SetFromDataElement(second_control_point.GetDataElement(at_weight.GetTag()));
		const double meterset_weight_per_spot = at_weight.GetValue() / double(NperSpot * n_spots_in_aperture);

#pragma omp parallel for
		for (int i_spot_xy = 0; (size_t)i_spot_xy < n_spots_in_aperture; i_spot_xy++)
		{
			const cl_float cur_thicc = thicc[i_spot_xy];
			const unsigned int i_group = n_total_spots + (i_spot_xy * j_W * NperSpot);

			// The random engine is NOT the source of our random problems!! AGA 15/08/2017
			std::minstd_rand0 g1((unsigned)std::clock() + i_spot_xy);  // minstd_rand0 is a standard linear_congruential_engine
			const float g1_max = float(g1.max());
			const float rand_to_angle = (2.0f * M_PI) / g1_max;
			const float rand_to_radius = 0.1f * mill / g1_max;

			for (size_t j_energy_step = 0; j_energy_step < j_W; j_energy_step++){
				for (size_t k = 0; k < NperSpot; k++){
					const unsigned int i_local = i_group + (j_energy_step * NperSpot) + k;

					// Dir(ection) is the opposite (unit)vector of that pointing from iso->source
					// Theta is the angle from z->phi and phi is angle from x->y
					// theta is gantry if z is towards anterior and phi is couch if y is cranial and x is lateral right
					float* tmp_float3 = new float[4]; // new must be used to ensure 16 byte alignment
					_mm_store_ps(tmp_float3, direction); // x-lateral   = -sin(theta)*cos(phi)
					// dir[i_local].s[1] = direction[1]; // y-AP        = -cos(theta)
					// dir[i_local].s[2] = direction[2]; // z-CC        = -sin(theta)*sin(phi)
					dir[i_local] = { tmp_float3[0], tmp_float3[1], tmp_float3[2] };

					const float rand_angle = float(g1()) * rand_to_angle;
					const float rand_radius = float(g1()) * rand_to_radius;

					const cl_float2 point = {
						spot_map[i_spot_xy].x + cos(rand_angle) * rand_radius,
						spot_map[i_spot_xy].y + sin(rand_angle) * rand_radius
					};

					const __m128 trans_point = { //   y=0, so only z of point
						A(0, 0) * point.x + A(0, 2) * point.y,
						A(1, 0) * point.x + A(1, 2) * point.y,
						A(2, 0) * point.x + A(2, 2) * point.y,
						0.0f // just to be 4 floats wide
					};

					_mm_store_ps(tmp_float3, _mm_add_ps(trans_point, beam_offset));
					pos[i_local] = { tmp_float3[0], tmp_float3[1], tmp_float3[2] };

					const cl_float de_comp = -bethe(nom_beam_energy + mod_energy_diff[j_energy_step]) * dedx_rel * cur_thicc;

					T[i_local] = nom_beam_energy + mod_energy_diff[j_energy_step] + de_comp; // mod_energy_diff already negated
					weight[i_local] = meterset_weight_per_spot * mod_energy_weigth[j_energy_step];
					// printf("E: %.3e, W: %.3e\n", T[i_local], weight[i_local]);
					delete[] tmp_float3;
				}
			}
		} // end for and omp parallel for

		delete[] mod_energy_diff;
		delete[] mod_energy_weigth;
		delete[] spot_map;

		n_total_spots += n_spots_in_aperture * j_W * NperSpot; // for enabeling full parallelism
	}

	return n_total_spots;
}

template<size_t N_per_spot>
size_t initSourceFromDicom(std::vector<std::string> plan_paths, cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight){
	size_t total_spots = 0;

	for (size_t i = 0; i < plan_paths.size(); i++){
		const char * dicom_path = plan_paths[i].c_str(); // args_info.plan_arg[i];
		// check dicom integrity
		gdcm::Reader reader;
		reader.SetFileName(dicom_path);
		if (!reader.Read())
		{
			std::cout << "Reading dicom plan failed!" << std::endl;
			return 0;
		}
		gdcm::File &file = reader.GetFile();
		gdcm::DataSet &ds = file.GetDataSet();
		//if (ds.FindDataElement(gdcm::Tag(0x10, 0x20)))
		//	const gdcm::DataElement &pt_id = ds.GetDataElement(gdcm::Tag(0x10, 0x20));

		const gdcm::DataElement &beam_seq_tag = ds.GetDataElement(gdcm::Tag(0x300a, 0x3a2));
		gdcm::SmartPointer<gdcm::SequenceOfItems> beam_seq = beam_seq_tag.GetValueAsSQ();

		const gdcm::Item &it_beams = beam_seq->GetItem(1);
		gdcm::Attribute< 0x300a, 0x308> at_scan_mode;
		at_scan_mode.SetFromDataElement(it_beams.GetDataElement(at_scan_mode.GetTag()));
		const char* scan_mode = at_scan_mode.GetValue();

		if (!strcmp(scan_mode, "NONE")){ // passive scatter
			size_t total_spots_cur = initFromPassiveScatter<N_per_spot>(beam_seq, T, pos, dir, weight, total_spots);
			std::cout << total_spots_cur << std::endl;
		}
		else // spot scanning
		{
			size_t total_spots_cur = initFromSpotScanning<N_per_spot>(beam_seq, T, pos, dir, weight, total_spots);
			std::cout << total_spots_cur << std::endl;
		}
		if (total_spots == 0)
			std::cout << "Something went wrong!! No points returned from dicom reader!!" << std::endl;
	}
	return total_spots;
}

// A function to initialize source protons. Should be replaced by real beams.
void initSource(cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight){
	int seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 g1((unsigned)seed1);  // minstd_rand0 is a standard linear_congruential_engine
	std::fill_n(T, N, 120.0f);
	std::fill_n(weight, N, 1.0f);
	for (int i = 0; i < N; i++){
		pos[i].s[0] = 5 * float(g1()) / g1.max(); // -15; // -15 to -10
		pos[i].s[1] = 120; // -20;   // ^--- between 0 and the largest possible max = 2147483646
		pos[i].s[2] = 5 * float(g1()) / g1.max(); // +25; //  25 to  30
	} //                                ^------- UP TO the largest possible max = 2147483646
	const cl_float3 temp2 = { 0.0f, 1.0f, 0.0f };
	std::fill_n(dir, N, temp2);
}

#ifdef GPMC_UI_H

QDir SaveUSHORTAsSHORT_DICOM_gdcmITK(itk::Image<unsigned short, 3>::Pointer& spImg, QString& strPatientID, QString& strPatientName, QString& strPathTargetDir)
{
	if (!spImg)
		return "";

	typedef itk::MinimumMaximumImageCalculator <UShortImageType> ImageCalculatorFilterType;
	ImageCalculatorFilterType::Pointer imageCalculatorFilter = ImageCalculatorFilterType::New();
	imageCalculatorFilter->SetImage(spImg);
	imageCalculatorFilter->Compute();

	double minVal = (double)(imageCalculatorFilter->GetMinimum());
	double maxVal = (double)(imageCalculatorFilter->GetMaximum());

	short outputMinVal = (short)(minVal - 1024);
	short outputMaxVal = (short)(maxVal - 1024);

	std::cout << "Output Min and Max Values are	" << outputMinVal << "	" << outputMaxVal << std::endl;

	typedef itk::RescaleIntensityImageFilter<UShortImageType, ShortImageType> RescaleFilterType;
	RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
	spRescaleFilter->SetInput(spImg);
	spRescaleFilter->SetOutputMinimum(outputMinVal);
	spRescaleFilter->SetOutputMaximum(outputMaxVal);
	spRescaleFilter->Update();

	ShortImageType::Pointer spShortImg = spRescaleFilter->GetOutput();

	QDir newDirPath(strPathTargetDir + "/" + strPatientID + strPatientName + "_DCM");

	QDir dirNew(newDirPath);
	if (!newDirPath.exists()) {
		newDirPath.mkdir(".");
	}
	else {
		if (newDirPath.removeRecursively()) {
			QDir dirReNew(newDirPath);
			dirReNew.mkdir(".");
		}
	}
	typedef itk::Image<unsigned short, 2> OutputImageType; //because dicom is one 2d image for each slice-file
	typedef itk::GDCMImageIO                ImageIOType;
	typedef itk::NumericSeriesFileNames     NamesGeneratorType;

	UShortImageType::RegionType region = spShortImg->GetLargestPossibleRegion();
	UShortImageType::IndexType  start = region.GetIndex();
	UShortImageType::SizeType   size = region.GetSize();

	ImageIOType::Pointer gdcmIO = ImageIOType::New();
	itk::MetaDataDictionary & dict = gdcmIO->GetMetaDataDictionary();
	std::string value;
	value = "CT";
	itk::EncapsulateMetaData<std::string>(dict, "0008|0060", value); // Modality
	value = "DERIVED\\SECONDARY\\AXIAL"; // This is virtually always correct when using ITK to write an image
	itk::EncapsulateMetaData<std::string>(dict, "0008|0008", value); // Image Type
	value = "SI";
	itk::EncapsulateMetaData<std::string>(dict, "0008|0064", value); // Conversion Type
	double value_double = spShortImg->GetSpacing()[2];
	std::ostringstream strs;
	strs << value_double;
	value = strs.str();
	// std::cout << "slice spacing: " + value << std::endl;
	itk::EncapsulateMetaData<std::string>(dict, "0018|0050", value); // SliceThickness
	itk::EncapsulateMetaData<std::string>(dict, "0018|0088", '-' + value); // SpacingBetweenSlices

	gdcm::UIDGenerator stduid;
	std::string studyUID = stduid.Generate();
	// std::cout << studyUID << std::endl;
	itk::EncapsulateMetaData<std::string>(dict, "0020|000d", studyUID);

	NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
	namesGenerator->SetStartIndex((unsigned)start[2]);
	namesGenerator->SetEndIndex(start[2] + size[2] - 1);
	namesGenerator->SetIncrementIndex(1);
	namesGenerator->SetSeriesFormat(newDirPath.absolutePath().toStdString() + "/CT." + studyUID + ".%d.dcm");

	typedef itk::ImageSeriesWriter<ShortImageType, OutputImageType> SeriesWriterType;
	SeriesWriterType::Pointer seriesWriter = SeriesWriterType::New();
	seriesWriter->SetInput(spShortImg);
	seriesWriter->SetImageIO(gdcmIO);
	seriesWriter->SetFileNames(namesGenerator->GetFileNames());

	try
	{
		seriesWriter->Update();
	}
	catch (itk::ExceptionObject & excp)
	{
		std::cerr << "Exception thrown while writing the series " << std::endl;
		std::cerr << excp << std::endl;
		return "";
	}
	// std::cout << "Alledgedly writing the series was successful to dir: " << newDirPath.toStdString() << std::endl;
	return newDirPath;
}

#endif