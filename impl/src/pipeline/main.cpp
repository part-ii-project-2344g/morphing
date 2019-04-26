#include"main.h"
#include"../common/typedefs.h"
#include"../common/constants.h"
#include"../embedding/io.h"
#include"../embedding/geometry_utils.h"
#include"../overlay/main.h"
#include"../interpolation/main.h"
#include"../measurement/main.h"
#include"../measurement_v/main.h"

#include <iostream>
#include <fstream>


/////////////////////////////////////////////////////////////////////////

#include"../embedding/main.h"
#include"../overlay/main.h"
#include"../interpolation/main.h"
#include"../measurement/main.h"
#include"../feature_alignment/main.h"
#include"../auto_features/main.h"
#include"../direct_optimization/main.h"

#include"../embedding/io.h"
#include"../common/constants.h"

#include"../auto_features/simple_features.h" // testing

#include <cstring>
#include <iostream>

#include <CGAL/Timer.h>

#include <windows.h>
#include <stdio.h>

using namespace std;

bool file_exists(const std::string& name) {
	ifstream f(name.c_str());
	return f.good();
}

int main_master(int argc, char **argv) {

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
		else {
			std::cout << "Usage: <subprogram-name> <args-for-subprogram>.\nAvailable subprograms:\nembedding\noverlay\ninterpolation\nmeasurement\nfeature_alignment\nauto_features\ndirect_optimization\n\n";
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
				&arg6[0u] };
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
			std::cout << "Usage: <subprogram-name> <args-for-subprogram>.\nAvailable subprograms:\nembedding\noverlay\ninterpolation\nmeasurement\npipeline\nfeature_alignment\nauto_features\n" << std::endl;
			ask_before_termination();
			return -1;
		}
	}
	subprogram_timer.stop();
	std::cout << "\n\n------------------------------------------\n[MASTER] Running time of the '" << subprogram_name << "' subprogram: " << subprogram_timer.time() << "s." << std::endl;
	return 0;
}


/////////////////////////////////////////////////////////////////////////









#define TEST true

#if TEST
int main_pipeline(int argc, char **argv) {
	CGAL::Verbose_ostream vout(true);
	vout << "\nVerbosity on." << endl;
	Polyhedron mesh_1, mesh_2;
	std::vector<Point> points_1{};
	std::vector<Point> points_2{};
	std::vector<std::vector<std::size_t>> faces_1{};
	std::vector<std::vector<std::size_t>> faces_2{};
	std::string input_path_1 = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\banana3_lowres.obj";
	std::string input_path_2 = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Deer\\fixed_deer3_lowres.obj";

	Polyhedron emesh_1, emesh_2;
	std::vector<Point> epoints_1{};
	std::vector<Point> epoints_2{};
	std::vector<std::vector<std::size_t>> efaces_1{};
	std::vector<std::vector<std::size_t>> efaces_2{};
	std::string einput_path_1 = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_banana3_lowres[embedding].obj";
	std::string einput_path_2 = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Deer\\_fixed_deer3_lowres[embedding].obj";

	std::string etemp_path_1 = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_banana3_lowres[rot_temp][embedding].obj";
	std::string etemp_path_1_read = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_banana3_lowres[rot_temp][embedding].obj";
	std::string temp_path_1 = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_banana3_lowres[rot_temp].obj";
	std::string temp_path_1_read = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_banana3_lowres[rot_temp].obj";
	std::string temp_overlay_path = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_banana3_lowres[rot_temp][overlay].obj";
	std::string temp_overlay_path_impersonation1 = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_banana3_lowres[rot_temp][overlay][impers_1].obj";
	std::string temp_overlay_path_impersonation2 = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_banana3_lowres[rot_temp][overlay][impers_2].obj";
	std::string temp_base_path = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models\\Banana\\_banana3_lowres[rot_temp][base_path].obj";

	std::string frames = "30";


	// Read meshes and embeddings.
	if (!input_mesh(points_1, faces_1, mesh_1, vout, input_path_1) || !input_mesh(points_2, faces_2, mesh_2, vout, input_path_2)
		|| !input_mesh(epoints_1, efaces_1, emesh_1, vout, einput_path_1) || !input_mesh(epoints_2, efaces_2, emesh_2, vout, einput_path_2)) {
		cerr << "Error reading meshes.\n";
		ask_before_termination();
		return -1;
	}
	std::string output_path = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\pipeline_measurements\\banana3_lowres@fixed_deer3_lowres[2]\\data.csv";

	// Prepare out file.
	ostream* p_out = &cout;
	ofstream out;
	out.open(output_path);
	p_out = &out;
	if (!*p_out) {
		cerr << "ERROR: cannot open file for writing: " << output_path << endl;
		ask_before_termination();
		return -1;
	}

	// Write the header.
	*p_out << "z,y,M,V" << endl;


	// Gather and write datapoints.
	int steps_x = 12;
	int steps_y = 12;

	// Count total_steps for debug purposes.
	int debug_total_steps = 0;
	int debug_steps_done = 0;
	for (int xi = -180; xi <= 180; xi += 2 * 180 / steps_x) {
		for (int yi = -180; yi <= 180; yi += 2 * 180 / steps_y) {
			debug_total_steps++;
		}
	}
	for (int zi = -180; zi <= 180; zi += 2*180 / steps_x) {
		for (int yi = -180; yi <= 180; yi += 2*180 / steps_y) {
			debug_steps_done++;
			double z = zi * CGAL_PI / 180.0; // convert to radians
			double y = yi * CGAL_PI / 180.0; // convert to radians
			Transformation t = rotZ(z) * rotY(y);

			Polyhedron emesh_1r;
			make_mesh_from_points_and_faces(epoints_1, efaces_1, emesh_1r);
			for (auto vertex_it = emesh_1r.vertices_begin(); vertex_it != emesh_1r.vertices_end(); vertex_it++) {
				(*vertex_it).point() = t.transform(vertex_it->point());
			}
			Polyhedron mesh_1r;
			make_mesh_from_points_and_faces(points_1, faces_1, mesh_1r);
			for (auto vertex_it = mesh_1r.vertices_begin(); vertex_it != mesh_1r.vertices_end(); vertex_it++) {
				(*vertex_it).point() = t.transform(vertex_it->point());
			}

			output_mesh(temp_path_1, mesh_1r, vout, false);
			output_mesh(etemp_path_1, emesh_1r, vout, false);

			// Perform the overlay
			char* args_overlay[7]{argv[0],
				"-v", 
				&temp_path_1_read[0u], // convert to char*
				&etemp_path_1_read[0u],
				&input_path_2[0u],
				&einput_path_2[0u],
				&temp_overlay_path[0u] };
			vout << "Running subprogram: overlay..." << endl << endl;
			main_overlay(7, args_overlay);

			// Perform the interpolation
			char* args_interpolation[6]{ argv[0],
				"-v",
				&frames[0u], // convert to char*
				&temp_overlay_path_impersonation1[0u],
				&temp_overlay_path_impersonation2[0u],
				&temp_base_path[0u] };
			vout << "Running subprogram: interpolation..." << endl << endl;
			main_interpolation(6, args_interpolation);

			// Perform the measurement(M)
			char* args_measurement[4]{ argv[0],
				"-v",
				&frames[0u], // convert to char*
				&temp_base_path[0u], };
			vout << "Running subprogram: measurement..." << endl << endl;
			double res_m = main_measurement_double(4, args_measurement);

			// Perform the measurement(V)
			char* args_measurement_v[4]{ argv[0],
				"-v",
				&frames[0u], // convert to char*
				&temp_base_path[0u], };
			vout << "Running subprogram: measurement_v..." << endl << endl;
			double res_v = main_measurement_v_double(4, args_measurement_v);

			// Write the datapoint.
			*p_out << zi << "," << yi << "," << res_m << "," << res_v << endl;

			// Debug.
			vout << "Written out the datapoint for (x=" << zi << ", y=" << yi << "). [" << (debug_steps_done*100)/debug_total_steps << "% done]" << endl;
		}
	}
	ask_before_termination();
	return 0;
}
#else
const int _roty_degrees = 80;
const std::string _roty80 = "_roty" + std::to_string(_roty_degrees);
std::string rotated_path(std::string path) {
	if (path.find_last_of("[") == std::string::npos) {
		return path.substr(0, path.find_last_of(".")) + _roty80 + path.substr(path.find_last_of("."));
	}
	else {
		return path.substr(0, path.find_last_of("[")) + _roty80 + path.substr(path.find_last_of("["));
	}
}
std::string rotated_name(std::string name) {
	return name + _roty80;
}
bool rotate_mesh_and_emb(std::string mesh_path, std::string emb_path) {
	std::cout << "Rotating mesh '" << mesh_path << "' and embedding '" << emb_path << "'..." << std::endl;
	Transformation rot = rotY(CGAL_PI * (double)_roty_degrees / 180.0);
	Polyhedron p, p2;
	std::vector<Point> points, points2;
	std::vector<std::vector<size_t>> faces, faces2;
	input_mesh(points, faces, p, CGAL::Verbose_ostream{}, mesh_path);
	input_mesh(points2, faces2, p2, CGAL::Verbose_ostream{}, emb_path);
	for (auto v = p.vertices_begin(); v != p.vertices_end(); v++) {
		v->point() = v->point().transform(rot);
	}
	for (auto v = p2.vertices_begin(); v != p2.vertices_end(); v++) {
		v->point() = v->point().transform(rot);
	}

	std::string out_mesh_path = rotated_path(mesh_path);
	std::string out_emb_path = rotated_path(emb_path);
	output_mesh(out_mesh_path, p, CGAL::Verbose_ostream{});
	output_mesh(out_emb_path, p2, CGAL::Verbose_ostream{});
	return true;
}
int main_pipeline(int argc, char **argv) {
	
	std::vector<std::tuple<std::string, std::string>> meshes{
		std::make_tuple("SimpleShapes", "simple_cube"),
		std::make_tuple("SimpleShapes", "simple_sphere"),
		std::make_tuple("SimpleShapes", "simple_tetrahedron_1"),
		std::make_tuple("SimpleShapes", "simple_tetrahedron_2"),
		std::make_tuple("Apple", "apple2_lowres"),
		std::make_tuple("Banana", "banana3_lowres"),
		std::make_tuple("Deer", "fixed_deer3_lowres"),
		std::make_tuple("Female", "female_lowres"),
		std::make_tuple("Fish", "fish2_lowres"),
		std::make_tuple("Horse", "horse2_lowres"),
		std::make_tuple("Male", "male_lowres"),
		std::make_tuple("Orange", "orange2_lowres"),
		std::make_tuple("Pear", "pear2_lowres"),
		std::make_tuple("Shark", "shark2_lowres"),
		std::make_tuple("Stones", "stone_long_lowres"),
		std::make_tuple("Stones", "stone_round_lowres")
	};

	size_t N = meshes.size();
	for (size_t i = 0; i < N-1; i++) {
		for (size_t j = i + 1; j < N; j++) {
			std::string folder1 = std::get<0>(meshes[i]);
			std::string folder2 = std::get<0>(meshes[j]);
			std::string filename1 = std::get<1>(meshes[i]);
			std::string filename2 = std::get<1>(meshes[j]);

			std::string models_folder = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\models";
			std::string mesh_1_path = models_folder + "\\" + folder1 + "\\"  + filename1 + ".obj";
			std::string emb_1_path  = models_folder + "\\" + folder1 + "\\_" + filename1 + "[embedding].obj";
			std::string mesh_2_path = models_folder + "\\" + folder2 + "\\"  + filename2 + ".obj";
			std::string emb_2_path  = models_folder + "\\" + folder2 + "\\_" + filename2 + "[embedding].obj";

			// Generate rotated meshes and embeddings, if not done yet.
			if (!(file_exists(rotated_path(mesh_1_path))) || !(file_exists(rotated_path(emb_1_path))))
				rotate_mesh_and_emb(mesh_1_path, emb_1_path);
			if (!(file_exists(rotated_path(mesh_2_path))) || !(file_exists(rotated_path(emb_2_path))))
				rotate_mesh_and_emb(mesh_2_path, emb_2_path);

			for (auto filenames_to_use : std::vector<std::pair<std::string, std::string>>{
				std::make_pair(filename1, filename2),
				std::make_pair(rotated_name(filename1), filename2)
				}) 
			{
				std::string filename1 = filenames_to_use.first;  //redefinition
				std::string filename2 = filenames_to_use.second; //redefinition

				char* args_none_1[16]{ argv[0],
					"-c",
					"feature_alignment",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"none",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_none_2[16]{ argv[0],
					"-c",
					"overlay",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"none",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_none_3[16]{ argv[0],
					"-c",
					"interpolation",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"none",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_none_4[16]{ argv[0],
					"-c",
					"measurement",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"none",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_none_5[16]{ argv[0],
					"-c",
					"measurement_v",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"none",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };

				char* args_poisson_0[16]{ argv[0],
					"-c",
					"auto_features",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"poisson",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_poisson_1[16]{ argv[0],
					"-c",
					"feature_alignment",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"poisson",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_poisson_2[16]{ argv[0],
					"-c",
					"overlay",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"poisson",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_poisson_3[16]{ argv[0],
					"-c",
					"interpolation",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"poisson",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_poisson_4[16]{ argv[0],
					"-c",
					"measurement",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"poisson",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_poisson_5[16]{ argv[0],
					"-c",
					"measurement_v",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"poisson",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };

				char* args_direct_1[16]{ argv[0],
					"-c",
					"direct_optimization",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"N/A",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_direct_2[16]{ argv[0],
					"-c",
					"overlay",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"direct",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_direct_3[16]{ argv[0],
					"-c",
					"interpolation",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"direct",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_direct_4[16]{ argv[0],
					"-c",
					"measurement",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"direct",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };
				char* args_direct_5[16]{ argv[0],
					"-c",
					"measurement_v",
					&folder1[0u], // convert to char*
					&folder2[0u],
					&filename1[0u],
					&filename2[0u],
					&filename1[0u],
					&filename2[0u],
					"direct",
					"40", "0.4", "3", "0.3", "1.0", "0.3" };

				std::string output_folder;
				bool skip = false;

				output_folder = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\morphs\\" + filename1 + "@" + filename2 + "@none";
				// Create output_folder.
				skip = false;
				if (!(CreateDirectory(output_folder.c_str(), NULL))) {
					if (ERROR_ALREADY_EXISTS != GetLastError()) {
						std::cerr << "Failed to create directory: " << output_folder << std::endl;
						return -1;
					}
					else {
						skip = true; // Skip if the directory already existed. Means we've failed before.
					}
				}
				if (!file_exists(output_folder + "\\_" + filename1 + "@" + filename2 + "@none-aligned.msr") && !skip) {
					std::cout << "{PIPELINE} Running " << filename1 + "@" + filename2 + "@none..." << std::endl;
					main_master(16, args_none_1);
					main_master(16, args_none_2);
					main_master(16, args_none_3);
					main_master(16, args_none_4); // msr
					main_master(16, args_none_5); // msrv
				}
				else {
					std::cout << "{PIPELINE} Skipping " << filename1 + "@" + filename2 + "@none..." << std::endl;
				}

				output_folder = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\morphs\\" + filename1 + "@" + filename2 + "@direct";
				// Create output_folder.
				skip = false;
				if (!(CreateDirectory(output_folder.c_str(), NULL))) {
					if (ERROR_ALREADY_EXISTS != GetLastError()) {
						std::cerr << "Failed to create directory: " << output_folder << std::endl;
						return -1;
					}
					else {
						skip = true;
					}
				}
				if (!file_exists(output_folder + "\\_" + filename1 + "@" + filename2 + "@direct-aligned.msr") && !skip) {
					std::cout << "{PIPELINE} Running " << filename1 + "@" + filename2 + "@direct..." << std::endl;
					main_master(16, args_direct_1);
					main_master(16, args_direct_2);
					main_master(16, args_direct_3);
					main_master(16, args_direct_4); // msr
					main_master(16, args_direct_5); // msrv
				}
				else {
					std::cout << "{PIPELINE} Skipping " << filename1 + "@" + filename2 + "@direct..." << std::endl;
				}

				output_folder = "D:\\!K\\UNIVERSITY\\YEAR3\\Dissertation\\morphing\\morphs\\" + filename1 + "@" + filename2 + "@poisson";
				// Create output_folder.
				skip = false;
				if (!(CreateDirectory(output_folder.c_str(), NULL))) {
					if (ERROR_ALREADY_EXISTS != GetLastError()) {
						std::cerr << "Failed to create directory: " << output_folder << std::endl;
						return -1;
					}
					else {
						skip = true;
					}
				}
				if (!file_exists(output_folder + "\\_" + filename1 + "@" + filename2 + "@poisson-aligned[embedding_1].obj") && !skip) { // Do not redo the ones that failed after alignment but before producing measurement.
					std::cout << "{PIPELINE} Running " << filename1 + "@" + filename2 + "@poisson..." << std::endl; 
					main_master(16, args_poisson_0);
					main_master(16, args_poisson_1);
					main_master(16, args_poisson_2);
					main_master(16, args_poisson_3);
					main_master(16, args_poisson_4); // msr
					main_master(16, args_poisson_5); // msrv
				}
				else {
					std::cout << "{PIPELINE} Skipping " << filename1 + "@" + filename2 + "@poisson..." << std::endl;
				}
			}
		}
	}

	return -1;
}
#endif
