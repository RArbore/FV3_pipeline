// Standard library
#include <iostream>
#include <string>
#include <vector>

// Unix
#include <sys/stat.h> 
#include <stdlib.h>
#include <unistd.h>

// MPI
#include <mpi.h>

// MODTRAN
#include <modtran/modlib/modlib.h>

// Constants
constexpr int MASTER    =  0; // Process with ID 0 is the master
constexpr int REQUEST   =  1; // Tag that requests are sent with over MPI
constexpr int REPLY     =  2; // Tag that replies are sent with over MPI
constexpr int DONE      = -1; // "Task id" sent to workers when there are no more tasks to do

// Forward declare the master routine
void master_routine(int num_proc, modtran::ModLib &modlib);

// Forward declare the worker routine
void worker_routine(int rank, modtran::ModLib &modlib);

// Check if a path is a directory
bool is_dir(const std::string& path) {
    struct stat stbuf;
    return (stat(path.c_str(), &stbuf) == 0 && S_ISDIR(stbuf.st_mode) != 0);
}

int main( int argc, char* argv[] ) {
    // Must provide 3 arguments: the input file, the output directory, and the data directory
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " input_file output_dir data_dir" << std::endl;
        return 0;
    }

    // Get current directory (will use later as MODTRAN appears to change directory we're in)
    char c_cwd[1024];
    getcwd(c_cwd, 1024);
    std::string cwd(c_cwd);
    cwd += "/";
    
    // Get some information from MPI (which process are we, how many processes are there)
    int rank, num_proc;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Convert the argument C-strings to std::strings
    std::string input_path(argv[1]);
    std::string output_path(argv[2]);
    std::string data_path(argv[3]);

    // Unique to this process, create a work directory
    std::string work_path;
    if (rank > 0) {
        work_path = output_path + "/" + std::to_string(rank) + std::string("/");
        if (mkdir(work_path.c_str(), 0777) == -1) {
            std::cerr << "ERROR: couldn't create temporary work directory " << work_path << std::endl;
            MPI_Finalize();
            return -1;
        }
    }
    else {
        work_path = output_path + "/";
    }
    
    // This script needs at least 2 processes (1 master, at least 1 worker)
    if (num_proc < 2) {
        std::cerr << "ERROR: at least two processes are required for execution." << std::endl;
        MPI_Finalize();
        return -1;
    }
    
    // Working directory must actually exist
    if (!is_dir(work_path)) {
        std::cerr << "ERROR: " << work_path << " does not exist or is not a directory." << std::endl;
        MPI_Finalize();
        return -1;
    }
    
    // Working directory must actually exist
    if (!is_dir(data_path)) {
        std::cerr << "ERROR: " << data_path << " does not exist or is not a directory." << std::endl;
        MPI_Finalize();
        return -1;
    }
    
    // Initialize the MODTRAN instance with the input file
    modtran::ModLib modlib;
    try {
        if (input_path.find( ".tp5" ) != std::string::npos ||
	    input_path.find( ".TP5" ) != std::string::npos)
	    modlib.inputTape5File(input_path.c_str());
        else
	    modlib.inputJsonFile(input_path.c_str());
    }
    catch (...) {
        std::cerr << "ERROR: a failure occurred while reading the input file " << input_path << std::endl
		  << (int) modlib.statusCode() << " " << modlib.statusMessage() << std::endl;
        MPI_Finalize();
        return -1;
    }
    
    // Attempt to load the working directory into MODTRAN
    if (!modlib.pathLocal(work_path.c_str())) {
        std::cerr << "ERROR: failed to set the working directory to " << work_path << std::endl;
        MPI_Finalize();
        return -1;
    }
    
    // Attempt to load the data directory into MODTRAN
    if (!modlib.pathData(data_path.c_str())) {
        std::cerr << "ERROR: failed to set the data directory to " << work_path << std::endl;
        MPI_Finalize();
        return -1;
    }
   
    // If we are the master, go into the master routine; if not, we go into the worker routine
    if (rank == MASTER) master_routine(num_proc, modlib);
    else worker_routine(rank, modlib);
    
    // Cleanup
    modlib.reset();
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    if (rank > 0) {
        std::string full_work_path = work_path[0] == '/' ? work_path : cwd + work_path;
        std::string full_output_path = output_path[0] == '/' ? output_path : cwd + output_path;
        std::string copy_cmd("mv " + full_work_path + "* " + full_output_path + "/");
        std::string remove_cmd("rm -r " + full_work_path);
        system(copy_cmd.c_str());
        system(remove_cmd.c_str());
    }
    return 0;
}

// Routine for master process
void master_routine(int num_proc, modtran::ModLib &modlib) {
    const int num_task = (int) modlib.caseCount();
    MPI_Status stat;
    int worker_id = 0, task_id = 0;

    // While there are tasks to complete, send task ids to available workers
    while (task_id < num_task) {
	MPI_Recv(&worker_id, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST, MPI_COMM_WORLD, &stat);
	MPI_Send(&task_id, 1, MPI_INT, worker_id, REPLY, MPI_COMM_WORLD);
	std::cout << "Executing case " << task_id << " on process " << worker_id << std::endl;
	++task_id;
    }

    // All tasks have been completed, so tell workers to stop (send DONE)
    for (int index = 0; index < num_proc-1; ++index) {
	MPI_Recv(&worker_id, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST, MPI_COMM_WORLD, &stat);
	MPI_Send(&DONE, 1, MPI_INT, worker_id, REPLY, MPI_COMM_WORLD);
    }
}

// Routine for worker processes
void worker_routine(int rank, modtran::ModLib &modlib) {
    MPI_Status stat;
    int task_id = 0;

    // Make myself available, get a task id
    MPI_Send(&rank, 1, MPI_INT, MASTER, REQUEST, MPI_COMM_WORLD);
    MPI_Recv(&task_id, 1, MPI_INT, MASTER, REPLY, MPI_COMM_WORLD, &stat);

    // While we haven't been told to finish, we keep reading task ids
    while (task_id != DONE) {
	ENModStatus case_stat = modlib.caseStatus(task_id);

	// If case still needs to be computed
	if (case_stat == STAT_INIT) {
	    // Modify some of the file output options of the case (add task id to name, add JSON output)
	    //ModInput* task_input = modlib.caseInput(task_id);
	    //std::string new_name = std::string(task_input->fileopt.flroot) + "_" + std::to_string(task_id);
	    //modtran::ModlibUtil::freeString(task_input->fileopt.flroot);
	    //task_input->fileopt.flroot = modtran::ModlibUtil::allocString(new_name.c_str());
	    //task_input->fileopt.jsonprnt = modtran::ModlibUtil::allocString(new_name.c_str());
	    //task_input->fileopt.jsonopt = WRT_ALL;
	    //task_input->fileopt.nofile = FC_NOFILES;

	    // Execute the case
            modlib.executeCase(task_id);
            modlib.caseOutputWrite(task_id);
	    modlib.caseOutputClear(task_id);
            if (modlib.statusCode() || modlib.caseStatus(task_id)) {
                std::cout << "Something went wrong, so here's the status message, status warning, case message, and case warning output from MODTRAN:\n"
                          << modlib.statusMessage() << std::endl
                          << modlib.statusWarning() << std::endl
                          << modlib.caseMessage(task_id) << std::endl
                          << modlib.caseWarning(task_id) << std::endl;
            }
	}

	// Make myself available, get a task id
	MPI_Send(&rank, 1, MPI_INT, MASTER, REQUEST, MPI_COMM_WORLD);
	MPI_Recv(&task_id, 1, MPI_INT, MASTER, REPLY, MPI_COMM_WORLD, &stat);
    }
}
