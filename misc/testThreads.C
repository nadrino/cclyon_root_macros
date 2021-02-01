


void *thread_fill_tree(void *i_thread_void_ptr);


void testThreads(){

  cout << "Starting testThreads..." << endl;

  vector<TThread*> tthreads(25); // 25 runs
  int nbAvailableSlots = 8-1; // 8 threads - 1 (this one)
  std::vector<int> jobStatusList(tthreads.size(), 0); // 0 -> not started, 1 -> started, 2 -> finished

  for(int i_run = 0 ; i_run < int(tthreads.size()) ; i_run++){

    string thread_name = "tree_filling_thread_" + to_string(i_run);

    while(nbAvailableSlots <= 0){
      std::this_thread::sleep_for(std::chrono::seconds(5));  // check every 5 sec
      for(int iThread = 0 ; iThread < int(jobStatusList.size()) ; iThread++){
        if(jobStatusList[iThread] == 1){     // if it is running
          if(tthreads[iThread]->GetState() >= TThread::kTerminatedState){ // if it's terminated but not Joined yet
            cout << "Thread" << iThread << " just finished, freeing one slot..." << endl;
            tthreads[iThread]->Join();       // wait for it to finish
            jobStatusList[iThread] = 2;      // tag it as finished
            nbAvailableSlots++;           // should leave the while loop and launch the next thread
            break;                        // leave this for loop
          }
          else{
            cout << "Thread" << iThread << " is still running: tthreads[iThread]->GetState() = " << tthreads[iThread]->GetState() << std::endl;
          }
        } // if thread running
      } // for
    } // while

    cout << "DISPATCHING NEW THREAD (" << i_run + 1 << "/" << tthreads.size() << endl;
    tthreads[i_run] = new TThread(thread_name.c_str(), thread_fill_tree, (void *) &i_run);
    tthreads[i_run]->Run();
    jobStatusList[i_run] = 1; // is running
    nbAvailableSlots--; // a slot is taken

    std::this_thread::sleep_for(std::chrono::seconds(2));

  }

  cout << "Ending testThreads..." << endl;

}


void *thread_fill_tree(void *i_thread_void_ptr){

  int i_thread = *((int *)i_thread_void_ptr);
  cout << "Starting thread #" << i_thread << " -> WAITING 20s" << endl;

  std::this_thread::sleep_for(std::chrono::seconds(20));

  cout << "Ending thread #" << i_thread << endl;

  return NULL;

}
