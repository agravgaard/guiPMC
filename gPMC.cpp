/*
 *
 * file  gPMC.cpp
 * brief program for using goPMC from commandline.
 * specifically made for use with CBCTrecon
 *
 * author Andreas Gravgaard Andersen
 *
 * last update on 17/8/2017
 *
 */
#include <algorithm>
#include <random>
#include <iostream>
#include <chrono>
#include <random>
#include <vector>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <itkEuler3DTransform.h>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include <itkImageSource.h>
#define NDOSECOUNTERS 1
#include "goPMC.h"

#define N_PER_SPOT 100

#include "range_modulator_data.hxx"
#include "gPMC_commandline.hxx"
#include "gPMC_dcm_tools.hxx"

// prototypes, real functions is defined after main.
void write_dose_to_mha(std::vector<cl_float> dose, args_info_gPMC &args_info);
void print_sum_and_mean(std::vector<cl_float> dose);
std::tuple<std::vector<cl_float>, std::vector<cl_float>> simulate(args_info_gPMC &args_info, const size_t N_dicom);

int main(int argc, char * argv[])
{
	GGO(gPMC, args_info);

	std::string lut_str(!args_info.lut_given ? "../lut" : args_info.lut_arg);
	std::string dicom_path(args_info.dir_arg);
	std::string plan_path(args_info.plan_arg);

	if (args_info.verbose_flag)
		std::cout << dicom_path << std::endl;

	// mcEngine.__autoclassinit2(651384); // just for hacks

	// Initialize source protons with arrays of energy (T), position (pos), direction (dir) and weight (weight) of each proton.
	// Position and direction should be defined in Dicom CT coordinate.

	std::cout << plan_path << "\nReading Dicom RT plan..." << std::endl;
	size_t N_dicom = GetNFromDicom<N_PER_SPOT>(plan_path);

	if (N_dicom == 69)
		return -1; // yes, 69 is a joke, but 69 is a highly improbable and specific number of spots.
	else if (args_info.verbose_flag)
		std::cout << " Dicom plan readable! " << N_dicom << " spots will be simulated." << std::endl;

	std::tuple<std::vector<cl_float>, std::vector<cl_float>> dose_tuple = simulate(args_info, N_dicom);

	std::vector<cl_float> doseMean = std::get<0>(dose_tuple);
	std::vector<cl_float> doseStd = std::get<1>(dose_tuple);

	for (size_t i = 1; i < (unsigned)args_info.batch_arg; i++){
		std::tuple<std::vector<cl_float>, std::vector<cl_float>> dose_tuple_tmp = simulate(args_info, N_dicom);
		std::vector<cl_float> doseMean_tmp = std::get<0>(dose_tuple_tmp);
		std::vector<cl_float> doseStd_tmp = std::get<1>(dose_tuple_tmp);

#pragma omp parallel for
		for (int j = 0; (size_t)j < doseMean.size(); j += 256){ // doseMean is always 256*256*N
			for (size_t j_local = (size_t)j; j_local < ((size_t)j + 256); j_local++)
				doseMean[j_local] += doseMean_tmp[j_local];
			// doseStd[j] += laplace(doseMean_tmp,j)^2 * doseStd_tmp[j]^2; (laplace(f, j) should be the derivative of f (a 3d-matrix) at j)
		}
		// sqrt(doseStd)
	}

	// Do something with doseMean and doseStd //
	if (args_info.verbose_flag){
		// Calculate mean of mean
		std::cout << "doseMean size: " << doseMean.size();
		print_sum_and_mean(doseMean);

		// Calculate mean of SD
		std::cout << "doseStd size: " << doseStd.size();
		print_sum_and_mean(doseStd);
	}
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	write_dose_to_mha(doseMean, args_info); // IO function

	// Clear the scoring counters in previous simulation runs.
	// mcEngine->clearCounter();
	// delete mcEngine;

	return 0;
}

bool check_sum(std::vector<cl_float> dose){
	cl_float sum = 0;

#pragma omp parallel for reduction(+:sum)
	for (int i = 0; (size_t)i < dose.size(); i += 256){ // doseMean is always 256*256*N
		cl_float local_sum = 0;

		for (size_t i_local = (size_t)i; i_local < ((size_t)i + 256); i_local++)
			local_sum += dose[i_local];

		sum += local_sum;
	}

	std::cout << "Sum: " << sum << std::endl;
	if (sum > 0.0f)
		return true;
	else
		return false;
}

std::tuple<cl::Platform, cl::Device> get_OpenCL_env(const char* hardware_arg){
	// Get OpenCL platform and device.
	cl::Platform platform;
	cl::Platform::get(&platform);

	std::vector<cl::Device> devs;

	if (!strcmp(hardware_arg, "gpu"))
		std::cout << "Getting device GPU returned: " << platform.getDevices(CL_DEVICE_TYPE_GPU, &devs) << std::endl;
	else if (!strcmp(hardware_arg, "cpu"))
		std::cout << "Getting devices CPU returned: " << platform.getDevices(CL_DEVICE_TYPE_CPU, &devs) << std::endl;
	else if (!strcmp(hardware_arg, "acc"))
		std::cout << "Getting devices ACCELERATOR returned: " << platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devs) << std::endl;
	else
		std::cout << "Getting devices DEFAULT returned: " << platform.getDevices(CL_DEVICE_TYPE_DEFAULT, &devs) << std::endl;

	cl::Device device;
	try{
		device = devs.at(0); // throws exception in contrary to []
		std::cout << "Using device: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "Error getting device: " << e.what() << std::endl;
		std::cerr << "Ususally means the program wasn't compiled for desired device!" << std::endl;
	}

	return std::make_tuple(platform, device);
}

goPMC::MCEngine* create_mcEngine(const char* hw_arg, const std::string dcm_dir, const std::string lut_str, const bool verbose){
	goPMC::MCEngine* mcEngine = new goPMC::MCEngine;

	std::tuple<cl::Platform, cl::Device> PlatDev_tuple = get_OpenCL_env(hw_arg);

	mcEngine->initializeComputation(std::get<0>(PlatDev_tuple), std::get<1>(PlatDev_tuple));
	// Read and process physics data.
	mcEngine->initializePhysics(lut_str); // Look Up Tables path, relative path accepted

	if (verbose)
		std::cout << "Physics initialised! Now reading dicom... " << std::endl;

	try{
		mcEngine->initializePhantom(dcm_dir); // "090737"); // "directoryToDicomData");
	}
	catch (const std::exception& e) {
		std::cerr << "Error reading dicom: " << e.what() << "\n"
			<< "Ususally means the implemented writer didn't give a compatible dicom image!" << std::endl;
		return mcEngine;
	}

	if (verbose)
		std::cout << "Dicom images read!" << std::endl;

	return mcEngine;
}

std::tuple<std::vector<cl_float>, std::vector<cl_float>>
call_gPMC(const size_t N_particles, const char* hardware_arg, const std::string dcm_dir, const std::string lut_dir, const bool verbose,
cl_float* T, cl_float3* pos, cl_float3* dir, cl_float* weight, std::string quantity)
{
	for (size_t i = 0; i < 10; i++) {
		std::vector<cl_float> doseMean, doseStd;
		goPMC::MCEngine* mcEngine = create_mcEngine(hardware_arg, dcm_dir, lut_dir, verbose); // Initialize OpenCL, physics and CT

		mcEngine->simulate(T, pos, dir, weight, N_particles, quantity);                       // Run simulation

		mcEngine->getResult(doseMean, doseStd);                                               // Get simulation results

		if ((check_sum(doseMean)) && (i > 0))
			std::cout << "Why did gPMC succeed on try number " << i << " and not the 1st?" << std::endl;

		mcEngine->clearCounter();                                                            // Clear the scoring counters in previous simulation runs. // No one knows
		delete mcEngine;

		if (check_sum(doseMean)) // check that there wasn't any odd side effects from mcEngine
			return std::make_tuple(doseMean, doseStd); // return if sum of doseMean > 0
	}
	std::vector<cl_float> errorVec = { 0.0f };
	return std::make_tuple(errorVec, errorVec);
}

std::tuple<std::vector<cl_float>, std::vector<cl_float>> simulate(args_info_gPMC &args_info, const size_t N_dicom)
{
	cl_float * T;
	cl_float3 * pos;
	cl_float3 * dir;
	cl_float * weight;

	if (!args_info.plan_given)
	{
		T = new cl_float[N];     //Energy(MeV?)= [120.0, ..., 120.0]
		pos = new cl_float3[N];  //Position    = [(5*rand_1-15, -20, 5*rand_1+25), ..., (5*rand_N-15, -20, 5*rand_N+25)]
		dir = new cl_float3[N];  //Direction   = [(0, 1, 0), ..., (0, 1, 0)] = y?          ^-------- 0 <= rand_X <= 1
		weight = new cl_float[N];//Weight      = [1.0, ..., 1.0]

		initSource(T, pos, dir, weight);
		if (args_info.verbose_flag)
			std::cout << "Source initialised, now simulating... ";
	}
	else
	{
		T = new cl_float[N_dicom];     //Energy(MeV?)= nominal control point energy
		pos = new cl_float3[N_dicom];  //Position    = Scan Spot Position transformed to gPMC coordinates
		dir = new cl_float3[N_dicom];  //Direction   = vector in unit sphere defined by the gantry and couch angle
		weight = new cl_float[N_dicom];//Weight      = Scan Spot Meterset Weights of control point

		std::vector<std::string> plan_paths;
		plan_paths.push_back(args_info.plan_arg);

		size_t N_result = initSourceFromDicom<N_PER_SPOT>(plan_paths, T, pos, dir, weight);
		if (N_result != N_dicom){
			std::cout << "\a" << "ONE OF THE COUNTERS ARE WRONG!!" << std::endl;
			std::cout << "dicom: " << N_dicom << ", result: " << N_result << std::endl;
		}

#pragma omp parallel for // (ignore signed/unsigned warning, omp only works with int)
		for (int i = 0; (size_t)i < N_dicom; i += N_PER_SPOT) // We know N_dicom is divisible by N_PER_SPOT by definition
		{
			for (size_t i_local = (size_t)i; i_local < ((size_t)i + N_PER_SPOT); i_local++){
				pos[i_local].s[0] *= 0.1f; // gPMC works in cm
				pos[i_local].s[1] *= 0.1f;
				pos[i_local].s[2] *= 0.1f;
			}
		}
	}

	// Choose a physics quantity to score for this simulation run.
	// Scoring quantity could be one of {DOSE2MEDIUM, DOSE2WATER, FLUENCE, LETD}.
	// LETD is dose weighted LET, to get dose averaged LET, divide it by DOSE2MEDIUM from another simulation run.
	//if (!strcmp(args_info.metric_arg, "dose2water"))
	std::string quantity("DOSE2WATER");
	//	quantity = "DOSE2WATER";
	if (!strcmp(args_info.metric_arg, "dose2medium"))
		quantity = "DOSE2MEDIUM";
	else if (!strcmp(args_info.metric_arg, "letd"))
		quantity = "LETD";
	else if (!strcmp(args_info.metric_arg, "fluence"))
		quantity = "FLUENCE";

	if (args_info.plan_given && args_info.verbose_flag){
		cl_float T_sum = 0.0f;
		cl_float W_sum = 0.0f;
		for (size_t idx = 0; idx < N_dicom; idx++){
			if (T[idx] < 0)
				std::cout << "T(" << idx << ") was less than zero: " << T[idx] << std::endl;
			else
				T_sum += T[idx];
			if (weight[idx] < 0)
				std::cout << "Weight(" << idx << ") was less than zero: " << weight[idx] << std::endl;
			else
				W_sum += weight[idx];
		}
		std::cout << "T sum: " << T_sum << " W sum: " << W_sum << std::endl;
	}
	size_t N_particles = (!args_info.plan_given ? N : N_dicom);
	std::string dcm_dir(args_info.dir_arg);
	std::string lut_dir(!args_info.lut_given ? "../lut" : args_info.lut_arg);
	// Run simulation.
	std::tuple<std::vector<cl_float>, std::vector<cl_float>> gPMC_output = call_gPMC(
		N_particles, args_info.hardware_arg, dcm_dir, lut_dir, args_info.verbose_flag, T, pos, dir, weight, quantity);

	// Get simulation results.
	std::vector<cl_float> doseMean = std::get<0>(gPMC_output);
	std::vector<cl_float> doseStd = std::get<1>(gPMC_output);

	if (args_info.verbose_flag){
		std::cout << "Dose Mean rigth after calc:" << std::endl;
		print_sum_and_mean(doseMean);
	}

	if (!strcmp(args_info.metric_arg, "letd")){
		quantity = "DOSE2MEDIUM";
		// Run simulation.
		std::tuple<std::vector<cl_float>, std::vector<cl_float>> gPMC_output = call_gPMC(
			N_particles, args_info.hardware_arg, dcm_dir, lut_dir, args_info.verbose_flag, T, pos, dir, weight, quantity);

		// Get simulation results.
		std::vector<cl_float> dmedMean = std::get<0>(gPMC_output);
		std::vector<cl_float> dmedStd = std::get<1>(gPMC_output);

#pragma omp parallel for
		for (int i = 0; (size_t)i < dmedMean.size(); i += 256){ // doseMean is always 256*256*N
			for (size_t i_local = (size_t)i; i_local < ((size_t)i + 256); i_local++)
				doseMean[i_local] /= dmedMean[i_local];
		}
	}
	delete[] T;
	delete[] pos;
	delete[] dir;
	delete[] weight;

	return std::make_tuple(doseMean, doseStd);
}

void print_sum_and_mean(std::vector<cl_float> dose){
	double sum = 0;

#pragma omp parallel for reduction(+:sum)
	for (int i = 0; (size_t)i < dose.size(); i += 256){ // doseMean is always 256*256*N
		double local_sum = 0;
		for (size_t i_local = (size_t)i; i_local < ((size_t)i + 256); i_local++)
			local_sum += dose[i_local];
		sum += local_sum;
	}

	printf(" sum: %.10e", sum);
	printf(" mean: %.10e\n", sum / dose.size());
}

void SetConstantImageSourceFromGgo(itk::Image<float, 3>::Pointer source, const args_info_gPMC &args_info)
{
	typedef itk::Image<float, 3> ImageType;
	const unsigned int Dimension = ImageType::GetImageDimension();

	ImageType::SizeType imageDimension;
	imageDimension.Fill((unsigned)args_info.dimension_arg[0]);
	if (args_info.dimension_given){
		imageDimension[0] = (unsigned)args_info.dimension_arg[0] / 2;
		imageDimension[1] = (unsigned)args_info.dimension_arg[1] / 2;
		imageDimension[2] = (unsigned)args_info.dimension_arg[2];
	}

	ImageType::SpacingType imageSpacing;
	imageSpacing.Fill(1.0);
	if (args_info.spacing_given){
		imageSpacing[0] = args_info.spacing_arg[0] * 2.0;
		imageSpacing[1] = args_info.spacing_arg[1] * 2.0;
		imageSpacing[2] = args_info.spacing_arg[2];
	}

	ImageType::PointType imageOrigin;
	for (unsigned int i = 0; i < Dimension; i++)
		imageOrigin[i] = imageSpacing[i] * (imageDimension[i] - 1) * -0.5;
	for (unsigned int i = 0; i < vnl_math_min(args_info.origin_given, Dimension); i++)
		imageOrigin[i] = args_info.origin_arg[i];

	ImageType::DirectionType imageDirection;
	if (args_info.direction_given)
		for (unsigned int i = 0; i < Dimension; i++)
			for (unsigned int j = 0; j < Dimension; j++)
				imageDirection[i][j] = args_info.direction_arg[i*Dimension + j];
	else
		imageDirection.SetIdentity();

	source->SetOrigin(imageOrigin);
	source->SetSpacing(imageSpacing);
	source->SetDirection(imageDirection);
	source->SetRegions(imageDimension);
	// source->SetConstant(0.);

	source->UpdateOutputInformation();
	source->Allocate();
}

void write_dose_to_mha(std::vector<cl_float> dose, args_info_gPMC &args_info){
	typedef itk::Image<float, 3> OutputImageType;
	OutputImageType::Pointer doseImage = OutputImageType::New();

	SetConstantImageSourceFromGgo(doseImage, args_info);

	//doseImage->Allocate();

	unsigned int i = 0;
	itk::ImageRegionIterator<OutputImageType> imIter(doseImage, doseImage->GetLargestPossibleRegion());
	while (!imIter.IsAtEnd())
	{
		imIter.Set(dose[i]);
		++imIter;
		++i;
	}

	if (args_info.verbose_flag)
		std::cout << "Writing output... " << std::endl;

	typedef  itk::ImageFileWriter<OutputImageType> WriterType;
	WriterType::Pointer outputWriter = WriterType::New();
	outputWriter->SetFileName(args_info.output_arg);
	outputWriter->SetInput(doseImage);
	outputWriter->Update();
}
