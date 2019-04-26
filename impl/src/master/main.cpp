#include"../embedding/main.h"
#include"../overlay/main.h"
#include"../interpolation/main.h"
#include"../measurement/main.h"
#include"../pipeline/main.h"
#include"../feature_alignment/main.h"
#include"../auto_features/main.h"
#include"../direct_optimization/main.h"
#include"../measurement_v/main.h"

#include"../embedding/io.h"
#include"../common/constants.h"

#include"../auto_features/simple_features.h" // testing

#include <cstring>
#include <iostream>

#include <CGAL/Timer.h>

// Example arguments:
/* 
embedding -v -s 1.097 0.11239 -0.06428 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Banana\banana3.obj"
embedding -v -s 1.097 0.11239 -0.06428 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Banana\banana3_lowres.obj"
embedding -v -s 1.72 0.0 0.0 0.0 1e-2 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Deer\fixed_deer3.obj"
embedding -v -s 1.72 0.0 0.0 0.0 1e-2 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Deer\fixed_deer3_lowres.obj"
embedding -v -s 1.52 0.34057 -0.0054 0.0 1e-5 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Horse\horse2.obj"
embedding -v -s 1.52 0.34057 -0.0054 0.0 1e-5 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Horse\horse2_lowres.obj"
embedding -v -s 1.87 0.9 0.0 -0.0 1e-5  "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Female\female.obj"
embedding -v -s 1.77 0.0 -0.0 0.0 1e-3 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Female\female_lowres.obj"
embedding -v -s 1.87 0.9 0.0 -0.0 1e-5  "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Male\male.obj"
embedding -v -s 1.77 0.0 -0.0 0.0 1e-3 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Male\male_lowres.obj"
embedding -v -s 1.4 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Fish\fish2.obj"
embedding -v -s 1.4 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Fish\fish2_lowres.obj"
embedding -v -s 2.5 0.305729 -0.036616 -0.006782 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Shark\shark2_lowres.obj"
embedding -v -s 1.0 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Orange\orange2.obj"
embedding -v -s 1.0 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Orange\orange2_lowres.obj"
embedding -v -s 1.0 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Apple\apple2_lowres.obj"
embedding -v -s 1.0 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Pear\pear2_lowres.obj"
embedding -v -s 1.0 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Pear\pear2.obj"
embedding -v -s 1.0 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Stones\stone_round.obj"
embedding -v -s 1.0 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Stones\stone_round_lowres.obj"
embedding -v -s 1.0 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Stones\stone_long.obj"
embedding -v -s 1.0 0.0 0.0 0.0 1e-4 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Stones\stone_long_lowres.obj"


==========> banana3 - fixed_deer3,   manual alignment: <==========
// Long form
feature_alignment -v "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual.feat" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Banana\banana3[embedding].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Deer\fixed_deer3[embedding].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned[embedding_1].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned[embedding_2].obj"
overlay -v "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Banana\banana3.obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned[embedding_1].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Deer\fixed_deer3.obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned[embedding_2].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\banana3@fixed_deer3@manual-aligned.obj"
interpolation -v 40 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned[impers_1].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned[impers_2].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned.obj"
measurement -v 40 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\frames\_banana3@fixed_deer3@manual-aligned.obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned.msr"

// Short form
-c feature_alignment Banana Deer banana3 fixed_deer3 banana3_lowres fixed_deer3_lowres poisson 40 0.4 3 0.3 1.0 0.3
-c overlay Banana Deer banana3 fixed_deer3 manual banana3_lowres fixed_deer3_lowres poisson 40 0.4 3 0.3 1.0 0.3
-c interpolation Banana Deer banana3 fixed_deer3 banana3_lowres fixed_deer3_lowres poisson 40 0.4 3 0.3 1.0 0.3
-c measurement Banana Deer banana3 fixed_deer3 banana3_lowres fixed_deer3_lowres poisson 40 0.4 3 0.3 1.0 0.3

// Directing stdout & stderr to output files:
overlay -v "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Banana\banana3.obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned[embedding_1].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\models\Deer\fixed_deer3.obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned[embedding_2].obj" "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned.obj" > "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\banana3@fixed_deer3@manual-aligned-overlay-log.mylog" 2>&1
measurement -v 40 "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\frames\_banana3@fixed_deer3@manual-aligned.obj" > "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana3@fixed_deer3@manual\_banana3@fixed_deer3@manual-aligned.msr" 2>&1

-c auto_features Banana Deer banana3 fixed_deer3 banana3_lowres fixed_deer3_lowres poisson 40 0.4 3 0.3 1.0 0.3
-c direct_optimization Banana Deer banana3_lowres fixed_deer3_lowres banana3_lowres fixed_deer3_lowres N/A 40 0.4 3 0.3 1.0 0.3

*/

// Entry point for the executable.
int main(int argc, char **argv) {

	// Tests.
	if (argc >= 2 && strcmp(argv[1], "test") == 0) {
		CGAL::Timer test_timer{}; test_timer.start();
		
		// Run test functions - they should fail via false assertions.
		test_path_underscore_functions();
		test_decimation();

		test_timer.stop(); std::cout << "All tests passed in " << test_timer.time() << "s." << std::endl;
		ask_before_termination();
		return 0;
	} // End of tests.

	// The core main method.
	CGAL::Timer subprogram_timer{};
	subprogram_timer.start();
	if (argc < 2) {
		std::cout << "Usage: <subprogram-name> <args-for-subprogram>.\nAvailable subprograms: \nembedding \noverlay \ninterpolation \nmeasurement \npipeline \n";
		ask_before_termination();
		return -1;
	}
	std::string subprogram_name = "error";
	if (strcmp(argv[1], "-c") != 0) {
		subprogram_name = argv[1];
		if (strcmp(argv[1], "embed") == 0 || strcmp(argv[1], "embedding") == 0) {
			for (int i = 1; i < argc - 1; i++)
				argv[i] = argv[i + 1];
			std::cout << "Running a subprogram: " << "embedding" << "..." << std::endl;
			main_embedding(argc - 1, argv);
		}
		else if (strcmp(argv[1], "overlay") == 0 || strcmp(argv[1], "merge") == 0) {
			for (int i = 1; i < argc - 1; i++)
				argv[i] = argv[i + 1];
			std::cout << "Running a subprogram: " << "overlay" << "..." << std::endl;
			main_overlay(argc - 1, argv);
		}
		else if (strcmp(argv[1], "interpolate") == 0 || strcmp(argv[1], "interpolation") == 0) {
			for (int i = 1; i < argc - 1; i++)
				argv[i] = argv[i + 1];
			std::cout << "Running a subprogram: " << "interpolation" << "..." << std::endl;
			main_interpolation(argc - 1, argv);
		}
		else if (strcmp(argv[1], "measurement") == 0 || strcmp(argv[1], "measure") == 0) {
			for (int i = 1; i < argc - 1; i++)
				argv[i] = argv[i + 1];
			std::cout << "Running a subprogram: " << "measurement" << "..." << std::endl;
			main_measurement(argc - 1, argv);
		}
		else if (strcmp(argv[1], "pipeline") == 0 || strcmp(argv[1], "pipeline") == 0) {
			for (int i = 1; i < argc - 1; i++)
				argv[i] = argv[i + 1];
			std::cout << "Running a subprogram: " << "pipeline" << "..." << std::endl;
			main_pipeline(argc - 1, argv);
		}
		else if (strcmp(argv[1], "feature_alignment") == 0 || strcmp(argv[1], "align") == 0) {
			for (int i = 1; i < argc - 1; i++)
				argv[i] = argv[i + 1];
			std::cout << "Running a subprogram: " << "feature_alignment" << "..." << std::endl;
			main_feature_alignment(argc - 1, argv);
		}
		else if (strcmp(argv[1], "auto_features") == 0 || strcmp(argv[1], "simple_features") == 0) {
			for (int i = 1; i < argc - 1; i++)
				argv[i] = argv[i + 1];
			std::cout << "Running a subprogram: " << "auto_features" << "..." << std::endl;
			main_simple_features(argc - 1, argv);
		}
		else if (strcmp(argv[1], "direct_optimization") == 0 || strcmp(argv[1], "direct") == 0) {
			for (int i = 1; i < argc - 1; i++)
				argv[i] = argv[i + 1];
			std::cout << "Running a subprogram: " << "direct_optimization" << "..." << std::endl;
			main_direct_optimization(argc - 1, argv);
		}
		else if (strcmp(argv[1], "measurement_v") == 0 || strcmp(argv[1], "measure_v") == 0) {
			for (int i = 1; i < argc - 1; i++)
				argv[i] = argv[i + 1];
			std::cout << "Running a subprogram: " << "measurement_v" << "..." << std::endl;
			main_measurement_v(argc - 1, argv);
		}
		else {
			std::cout << "Usage: <subprogram-name> <args-for-subprogram>.\nAvailable subprograms:\nembedding\noverlay\ninterpolation\nmeasurement\npipeline\nfeature_alignment\nauto_features\ndirect_optimization\nmeasurement_v\n\n";
			ask_before_termination();
			return -1;
		}
	}
	else if (strcmp(argv[1], "-c") == 0) { // Generate the right paths automatically.
		std::string root = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\";
		subprogram_name = argv[2]; // e.g. embedding
		std::string model_folder1 = argv[3]; // e.g. Banana
		std::string model_folder2 = argv[4]; // e.g. Deer
		std::string name1 = argv[5]; // e.g. banana3
		std::string name2 = argv[6]; // e.g. fixed_deer3
		std::string name1_lowres = argv[7]; // e.g. banana3_lowres
		std::string name2_lowres = argv[8]; // e.g. fixed_deer3_lowres
		std::string mode = argv[9]; // e.g. manual
		std::string frames = argv[10]; // e.g. 40
		std::string p_check = argv[11]; // e.g. 0.4
		std::string n_feat = argv[12]; // e.g. 3
		std::string r_p = argv[13]; // e.g. 0.3
		std::string dist_weight = argv[14]; // e.g. 1.0
		std::string r_f = argv[15]; // e.g. 0.3
		std::string combined = name1 + "@" + name2 + "@" + mode;
		std::string combined_direct_opt = name1 + "@" + name2 + "@" + "direct";

		if (subprogram_name == "overlay") {
			std::string arg2 = root + "models\\" + model_folder1 + "\\_" + name1 + ".obj";
			std::string arg3 = root + "morphs\\" + combined + "\\_" + combined + "-aligned[embedding_1].obj";
			std::string arg4 = root + "models\\" + model_folder2 + "\\_" + name2 + ".obj";
			std::string arg5 = root + "morphs\\" + combined + "\\_" + combined + "-aligned[embedding_2].obj";
			std::string arg6 = root + "morphs\\" + combined + "\\_" + combined + "-aligned.obj";

			char* args[7]{ argv[0],
				"-v",
				&arg2[0u], // convert to char*
				&arg3[0u],
				&arg4[0u],
				&arg5[0u],
				&arg6[0u] };

			std::cout << "Running a subprogram: " << "overlay" << "..." << std::endl;
			main_overlay(7, args);
		}
		else if (subprogram_name == "interpolate" || subprogram_name == "interpolation") {

			std::string arg2 = frames;
			std::string arg3 = root + "morphs\\" + combined + "\\_" + combined + "-aligned[impers_1].obj";
			std::string arg4 = root + "morphs\\" + combined + "\\_" + combined + "-aligned[impers_2].obj";
			std::string arg5 = root + "morphs\\" + combined + "\\_" + combined + "-aligned.obj";
			
			char* args[6]{ argv[0],
				"-v",
				&arg2[0u], // convert to char*
				&arg3[0u],
				&arg4[0u],
				&arg5[0u] };

			std::cout << "Running a subprogram: " << "interpolation" << "..." << std::endl;
			main_interpolation(6, args);
		}
		else if (subprogram_name == "measurement") {

			// measurement -v 40 
			// "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana2@fixed_deer2@manual\frames\_banana2@fixed_deer2@manual-aligned.obj" 
			// "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana2@fixed_deer2@manual\_banana2@fixed_deer2@manual-aligned.msr"

			std::string arg2 = frames;
			std::string arg3 = root + "morphs\\" + combined + "\\_" + combined + "-aligned.obj";
			std::string arg4 = root + "morphs\\" + combined + "\\_" + combined + "-aligned.msr";

			char* args[5]{ argv[0],
				"-v",
				&arg2[0u], // convert to char*
				&arg3[0u],
				&arg4[0u] };

			std::cout << "Running a subprogram: " << "measurement" << "..." << std::endl;
			main_measurement(5, args);
		}
		else if (subprogram_name == "measurement_v") {

			// measurement_v -v 40 
			// "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana2@fixed_deer2@manual\frames\_banana2@fixed_deer2@manual-aligned.obj" 
			// "D:\!K\UNIVERSITY\YEAR3\Dissertation\morphing\morphs\banana2@fixed_deer2@manual\_banana2@fixed_deer2@manual-aligned.msr"

			std::string arg2 = frames;
			std::string arg3 = root + "morphs\\" + combined + "\\_" + combined + "-aligned.obj";
			std::string arg4 = root + "morphs\\" + combined + "\\_" + combined + "-aligned.msrv";

			char* args[5]{ argv[0],
				"-v",
				&arg2[0u], // convert to char*
				&arg3[0u],
				&arg4[0u] };

			std::cout << "Running a subprogram: " << "measurement_v" << "..." << std::endl;
			main_measurement_v(5, args);
		}
		else if (subprogram_name == "feature_alignment") {

			// If there is no alignment, then just copy the embeddings to 'aligned embeddings' by passing empty string as feature path.
			std::string arg2 = (mode == "none") ? "" : (root + "morphs\\" + combined + "\\_" + combined + ".feat");
			std::string arg3 = root + "models\\" + model_folder1 + "\\_" + name1 + "[embedding].obj";
			std::string arg4 = root + "models\\" + model_folder2 + "\\_" + name2 + "[embedding].obj";
			std::string arg5 = root + "morphs\\" + combined + "\\_" + combined + "-aligned[embedding_1].obj";
			std::string arg6 = root + "morphs\\" + combined + "\\_" + combined + "-aligned[embedding_2].obj";

			char* args[7]{ argv[0],
				"-v",
				&arg2[0u], // convert to char*
				&arg3[0u],
				&arg4[0u],
				&arg5[0u],
				&arg6[0u]};
			std::cout << "Running a subprogram: " << "feature_alignment" << "..." << std::endl;
			main_feature_alignment(7, args);
		}
		else if (subprogram_name == "auto_features" || subprogram_name == "simple_features") {

			std::string arg2 = mode;
			std::string arg3 = root + "models\\" + model_folder1 + "\\_" + name1 + ".obj";
			std::string arg4 = root + "models\\" + model_folder2 + "\\_" + name2 + ".obj";
			std::string arg5 = root + "models\\" + model_folder1 + "\\_" + name1_lowres + ".obj";
			std::string arg6 = root + "models\\" + model_folder2 + "\\_" + name2_lowres + ".obj";
			std::string arg7 = root + "morphs\\" + combined + "\\_" + combined + ".feat";

			char* args[13]{ argv[0],
				"-v",
				&arg2[0u], // convert to char*
				&arg3[0u],
				&arg4[0u],
				&arg5[0u],
				&arg6[0u],
				&arg7[0u],
				&p_check[0u],
				&n_feat[0u],
				&r_p[0u],
				&dist_weight[0u],
				&r_f[0u] };

			std::cout << "Running a subprogram: " << "auto_features" << "..." << std::endl;
			main_simple_features(13, args);
		}
		else if (subprogram_name == "direct_optimization" || subprogram_name == "direct") {

			std::string arg2 = root + "models\\" + model_folder1 + "\\_" + name1 + ".obj";
			std::string arg3 = root + "models\\" + model_folder1 + "\\_" + name1 + "[embedding].obj";
			std::string arg4 = root + "models\\" + model_folder2 + "\\_" + name2 + ".obj";
			std::string arg5 = root + "models\\" + model_folder2 + "\\_" + name2 + "[embedding].obj";
			std::string arg6 = root + "morphs\\" + combined_direct_opt + "\\_" + combined_direct_opt + "[temp].obj"; // output base
			std::string arg7 = root + "morphs\\" + combined_direct_opt + "\\_" + combined_direct_opt + "-aligned[embedding_1].obj"; // output_1
			std::string arg8 = root + "morphs\\" + combined_direct_opt + "\\_" + combined_direct_opt + "-aligned[embedding_2].obj"; // output_2

			char* args[14]{ argv[0],
				"-v",
				&arg2[0u], // convert to char*
				&arg3[0u],
				&arg4[0u],
				&arg5[0u],
				&arg6[0u],
				&arg7[0u],
				&arg8[0u],
				&p_check[0u],
				&n_feat[0u],
				&r_p[0u],
				&dist_weight[0u],
				&r_f[0u] };

			std::cout << "Running a subprogram: " << "direct_optimization" << "..." << std::endl;
			main_direct_optimization(14, args);
		}
		else {
			std::cout << "Subprogram " << subprogram_name << " does not exist or does not support the '-c' option." << std::endl;
			std::cout << "Usage: <subprogram-name> <args-for-subprogram>.\nAvailable subprograms:\nembedding\noverlay\ninterpolation\nmeasurement\npipeline\nfeature_alignment\nauto_features\nmeasurement_v\n" << std::endl;
			ask_before_termination();
			return -1;
		}
	}
	subprogram_timer.stop();
	std::cout << "\n\n------------------------------------------\n[MASTER] Running time of the '" << subprogram_name << "' subprogram: " << subprogram_timer.time() << "s." << std::endl;
	return 0;
}








