{

  cout << "Using ROOT logon : $REPO_DIR/cclyon_root_macros/logon/rootlogon.C " << endl;

  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatFont(42);
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(0);

  gStyle->SetLabelFont(42,"xyz");
  gStyle->SetLabelSize(0.05,"xyz");
  gStyle->SetLabelOffset(0.015,"x");
  gStyle->SetLabelOffset(0.015,"y");

  gStyle->SetTitleFont(42);
  gStyle->SetTitleFontSize(0.06);
  gStyle->SetLegendFont(42);

  gStyle->SetTitleFont(42,"xyz");
  gStyle->SetTitleSize(0.06,"xyz");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleOffset(1.20,"y");

  gStyle->SetStripDecimals(kFALSE);

  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.07);
  gStyle->SetPadTopMargin(0.09);

  // gStyle->SetStatW(0.35);
  // gStyle->SetStatH(0.25);

  gStyle->SetPadTickX(kTRUE);
  gStyle->SetPadTickY(kTRUE);

  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);

  //  gStyle->SetPalette(1);

  gStyle->SetLineWidth(2);
  gStyle->SetHistLineWidth(3);
  gStyle->SetFuncWidth(3);
  gStyle->SetFrameLineWidth(2);

  gStyle->SetMarkerSize(1.2);

  if(gROOT->GetVersionInt() > 60000){
    std::cout << "Loading GenericToolbox namespace..." << std::endl;
    std::string REPO_DIR(getenv("REPO_DIR"));
    std::string s = ".include " + REPO_DIR + "/cclyon_root_macros/submodules/cpp-generic-toolbox/include";
    // R__ADD_INCLUDE_PATH((REPO_DIR + "/cclyon_root_macros/submodules/cpp-generic-toolbox/include").c_str());
    // gROOT->LoadMacro("/Users/ablanche/Documents/Work/Repositories/cclyon_root_macros/submodules/cpp-generic-toolbox/include/GenericToolbox.h");
    gROOT->ProcessLine(s.c_str());
    gROOT->ProcessLine("#include \"GenericToolbox.h\"");
    gROOT->ProcessLine("#include \"GenericToolboxRootExt.h\"");
  }

}


#include <sys/types.h>
#include <sys/stat.h>

namespace TToolBox {

  static string NORMAL = "\033[00m";
  static string ERROR    = "\033[1;31m<ERROR> \033[00m";
  static string INFO  = "\033[1;32m<INFO> \033[00m";
  static string WARNING   = "\033[1;33m<WARNING> \033[00m";
  static string ALERT = "\033[1;35m<ALERT> \033[00m";
  static TH1D* dummy_TH1D = nullptr;

  // Forward declarations
  std::vector<std::string> split_string(std::string input_string_, std::string delimiter_);
  std::string join_vector_string(std::vector<std::string> string_list_, std::string delimiter_, int begin_index_, int end_index_);

  ////////////////////////////
  // Src
  static double last_displayed_value = -1;
  static int verbosity_level = 0;
  static Color_t colorCells[10] = {kOrange+1, kGreen-3, kTeal+3, kAzure+7, kCyan-2, kBlue-7, kBlue+2, kOrange+9, kRed+2, kPink+9};

  void display_loading(int current_index_, int end_index_, string title_ = "", bool force_display_ = false) {

    int percent = int(round(double(current_index_) / end_index_ * 100.));
    if(force_display_ or current_index_ >= end_index_-1) {
      if(last_displayed_value != -1) cout << "\r" << title_ << " : " << 100 << "%" << endl;
      last_displayed_value = -1;
      return;
    }
    if(last_displayed_value == -1 or last_displayed_value < percent) {
      last_displayed_value = percent;
      cout << "\r" << title_ << " : " << percent << "%" << flush << "\r";
    }

  }
  void write_string_in_file(string file_path_, string string_to_write_){
    std::ofstream out(file_path_.c_str());
    out << string_to_write_;
    out.close();
  }

  void fix_TH2D_display(TH2D* histogram_){

    gPad->SetRightMargin(0.15);
    histogram_->GetZaxis()->SetTitleOffset(0.8);
    TPaletteAxis* pal = (TPaletteAxis*) histogram_->GetListOfFunctions()->FindObject("palette");
    // TPaletteAxis* pal = (TPaletteAxis*) histogram_->GetListOfFunctions()->At(0);
    if(pal != nullptr){
      pal->SetX1NDC(1 - 0.15 + 0.01);
      pal->SetX2NDC(1 - 0.15 + 0.05);
      pal->GetAxis()->SetMaxDigits(2);
      pal->Draw();
    }

  }
  void set_pallette_default(){
    gStyle->SetPalette(kBird);
  }
  void set_pallette_blue_red(){
    gStyle->SetPalette(kBlackBody);
    TColor::InvertPalette();
    // gStyle->SetPalette();
  }
  void set_pallette_orange(){
    gStyle->SetPalette(kDarkBodyRadiator);
  }

  bool do_path_is_valid(std::string path_){
    struct stat buffer{};
    return (stat (path_.c_str(), &buffer) == 0);
  }
  bool do_path_is_folder(std::string folder_path_){
    struct stat info{};
    stat( folder_path_.c_str(), &info );
    if( info.st_mode & S_IFDIR ) return true;
    else return false;
  }
  bool do_path_is_file(std::string file_path_) {
    if(do_path_is_valid(file_path_)){
      return not do_path_is_folder(file_path_);
    } else{
      return false;
    }
  }
  bool do_tfile_is_valid(TFile *input_tfile_, bool check_if_writable_=false){

      if(input_tfile_ == nullptr){
          if(verbosity_level >= 1) std::cerr << ERROR << "input_tfile_ is a nullptr" << std::endl;
          return false;
      }

      if(not input_tfile_->IsOpen()){
          if(verbosity_level >= 1) std::cerr << ERROR << "input_tfile_ = " << input_tfile_->GetName() << " is not opened." << std::endl;
          if(verbosity_level >= 1) std::cerr << ERROR << "input_tfile_->IsOpen() = " << input_tfile_->IsOpen() << std::endl;
          return false;
      }

      if(check_if_writable_ and not input_tfile_->IsWritable()){
          if(verbosity_level >= 1) std::cerr << ERROR << "input_tfile_ = " << input_tfile_->GetName() << " is not writable." << std::endl;
          if(verbosity_level >= 1) std::cerr << ERROR << "input_tfile_->IsWritable() = " << input_tfile_->IsWritable() << std::endl;
          return false;
      }

      return true;

  }
  bool do_tfile_is_valid(std::string input_file_path_){
    bool file_is_valid = false;
    if(do_path_is_file(input_file_path_)){
        auto old_verbosity = gErrorIgnoreLevel;
        gErrorIgnoreLevel = kFatal;
        auto* input_tfile = TFile::Open(input_file_path_.c_str(), "READ");
        if(do_tfile_is_valid(input_tfile)){
            file_is_valid = true;
            input_tfile->Close();
        }
        delete input_tfile;
        gErrorIgnoreLevel = old_verbosity;
    }
    return file_is_valid;
}
  bool do_string_in_vector(std::string str_, std::vector<std::string>& vector_){
    for(auto const &element : vector_){
      if(element == str_) return true;
    }
    return false;
  }
  bool do_string_contains_substring(std::string string_, std::string substring_){
    if(substring_.size() > string_.size()) return false;
    if(string_.find(substring_) != std::string::npos) return true;
    else return false;
  }
  bool do_string_starts_with_substring(std::string string_, std::string substring_){
    if(substring_.size() > string_.size()) return false;
    return (not string_.compare(0, substring_.size(), substring_));
  }
  bool do_string_ends_with_substring(std::string string_, std::string substring_){
    if(substring_.size() > string_.size()) return false;
    return (not string_.compare(string_.size() - substring_.size(), substring_.size(), substring_));
  }

  std::string get_folder_path_from_file_path(std::string file_path_){
    std::string folder_path;
    if(file_path_[0] == '/') folder_path += "/";
    auto splitted_path = split_string(file_path_, "/");
    folder_path += join_vector_string(
        splitted_path,
        "/",
        0,
        int(splitted_path.size()) - 1
        );
    return folder_path;
  }
  std::string get_filename_from_file_path(std::string file_path_){
    auto splitted_path = split_string(file_path_, "/");
    return splitted_path.back();
  }
  std::string get_filename_without_ext_from_file_path(std::string file_path_){
    auto file_name_elements = split_string(get_filename_from_file_path(file_path_), ".");
    return join_vector_string(file_name_elements, ".", 0, file_name_elements.size()-1);
  }
  std::string join_vector_string(std::vector<std::string> string_list_, std::string delimiter_, int begin_index_ = 0, int end_index_ = 0) {

    std::string joined_string;
    if(end_index_ == 0) end_index_ = int(string_list_.size());

    // circular permutation -> python style : tab[-1] = tab[tab.size - 1]
    if(end_index_ < 0 and int(string_list_.size()) > std::fabs(end_index_)) end_index_ = int(string_list_.size()) + end_index_;

    for(int i_list = begin_index_ ; i_list < end_index_ ; i_list++){
      if(not joined_string.empty()) joined_string += delimiter_;
      joined_string += string_list_[i_list];
    }

    return joined_string;
  }

  std::vector<double> get_log_binning(int n_bins_, double X_min_, double X_max_){
    vector<double> output(n_bins_+1); // add one extra bin for the boundary
    double xlogmin = TMath::Log10(X_min_);
    double xlogmax = TMath::Log10(X_max_);
    double dlogx   = (xlogmax-xlogmin)/((double)n_bins_);
    for( int i_bin = 0 ; i_bin <= n_bins_; i_bin++) {
      double xlog = xlogmin + i_bin*dlogx;
      output[i_bin] = TMath::Exp( TMath::Log(10) * xlog );
    }
    return output;
  }
  std::vector<double> get_linear_binning(int n_bins_, double X_min_, double X_max_){
    vector<double> output(n_bins_+1); // add one extra bin for the boundary
    double dx   = (X_max_-X_min_)/((double)n_bins_);
    for( int i_bin = 0 ; i_bin <= n_bins_; i_bin++) {
      double x = X_min_ + i_bin*dx;
      output[i_bin] = x;
    }
    return output;
  }

  std::vector<std::string> split_string(std::string input_string_, std::string delimiter_){

    std::vector<std::string> output_splited_string;

    const char *src = input_string_.c_str();
    const char *next = src;

    std::string out_string_piece = "";

    while ((next = std::strstr(src, delimiter_.c_str())) != NULL) {
      out_string_piece = "";
      while (src != next){
        out_string_piece += *src++;
      }
      output_splited_string.emplace_back(out_string_piece);
      /* Skip the delimiter_ */
      src += delimiter_.size();
    }

    /* Handle the last token */
    out_string_piece = "";
    while (*src != '\0')
      out_string_piece += *src++;

    output_splited_string.emplace_back(out_string_piece);

    return output_splited_string;

  }
  std::vector<std::string> get_list_of_entries_in_folder(std::string *folder_path_){

    if(not do_path_is_folder(*folder_path_)) return std::vector<std::string>();

    std::vector<std::string> entries_list;
    TSystemDirectory dir((*folder_path_).c_str(), (*folder_path_).c_str());
    TList *files = dir.GetListOfFiles();
    for(int i_entry = 0 ; i_entry < dir.GetListOfFiles()->GetEntries() ; i_entry++){
      string entry = dir.GetListOfFiles()->At(i_entry)->GetName();
      if(entry != "." and entry != ".."){
        entries_list.emplace_back(entry);
      }
    }
    return entries_list;

  }
  std::vector<std::string> get_list_of_subfolders_in_folder(std::string *folder_path_) {

    auto entries_list = get_list_of_entries_in_folder(folder_path_);
    std::vector<std::string> folders_list;
    for(int i_entry = 0 ; i_entry < int(entries_list.size()) ; i_entry++){
      if(do_path_is_folder(*folder_path_ + "/" + entries_list[i_entry]))
        folders_list.emplace_back(entries_list[i_entry]);
    }
    return folders_list;

  }
  std::vector<std::string> get_list_of_files_in_folder(std::string *folder_path_, std::string *files_extension_ = nullptr) {

    auto entries_list = get_list_of_entries_in_folder(folder_path_);
    std::vector<std::string> files_list;
    for(int i_entry = 0 ; i_entry < int(entries_list.size()) ; i_entry++){
      if(files_extension_ == nullptr or do_string_ends_with_substring(entries_list[i_entry], *files_extension_)){
        if(do_path_is_file(*folder_path_ + "/" + entries_list[i_entry]))
          files_list.emplace_back(entries_list[i_entry]);
      }
    }
    return files_list;

  }
  std::vector<std::string> get_list_of_files_in_subfolders(std::string *folder_path_, std::string *files_extension_ = nullptr){

    std::vector<std::string> output_file_paths;

    auto files_list = get_list_of_files_in_folder(folder_path_, files_extension_);
    for(int i_file = 0 ; i_file < int(files_list.size()) ; i_file++){
      output_file_paths.emplace_back(files_list[i_file]);
    }

    auto subfolders_list = get_list_of_subfolders_in_folder(folder_path_);
    for(int i_subfolder = 0 ; i_subfolder < int(subfolders_list.size()) ; i_subfolder++){
      string subfolde_full_path = *folder_path_ + "/" + subfolders_list[i_subfolder];
      auto subfiles_path = get_list_of_files_in_subfolders(&subfolde_full_path, files_extension_); // RECURSIVE
      for(int i_subfile = 0 ; i_subfile < int(subfiles_path.size()) ; i_subfile++){
         output_file_paths.emplace_back(subfolders_list[i_subfolder] + "/" + subfiles_path[i_subfile]);
      }
    }

    return output_file_paths;

  }
  std::vector<std::string> read_file(std::string file_path_){

    std::vector<std::string> res;
    std::string line;
    std::ifstream infile(file_path_.c_str());
    while (std::getline(infile, line))
    {
      res.emplace_back(line);
    }
    infile.close();

    return res;

  }

  // Overloaded (slower but easier to use)
  std::vector<std::string> get_list_of_entries_in_folder(std::string folder_path_){
    std::string folder_path = folder_path_;
    return get_list_of_entries_in_folder(&folder_path);
  }
  std::vector<std::string> get_list_of_subfolders_in_folder(std::string folder_path_){
    std::string folder_path = folder_path_;
    return get_list_of_subfolders_in_folder(&folder_path);
  }
  std::vector<std::string> get_list_of_files_in_folder(std::string folder_path_, std::string files_extension_ = ""){
    std::string *folder_path = &folder_path_;
    std::string *files_extension = nullptr;
    if(not files_extension_.empty()) files_extension = &files_extension_;
    return get_list_of_files_in_folder(folder_path, files_extension);
  }
  std::vector<std::string> get_list_of_files_in_subfolders(std::string folder_path_, std::string files_extension_ = ""){
    std::string *folder_path = &folder_path_;
    std::string *files_extension = nullptr;
    if(not files_extension_.empty()) files_extension = &files_extension_;
    return get_list_of_files_in_subfolders(folder_path, files_extension);
  }

  // Matrices/Vector Tools
  TH1D* get_TH1D_log_binning(string name_, string title_, int n_bins_, double X_min_, double X_max_){

    TH1D* output = nullptr;
    vector<double> xbins = get_log_binning(n_bins_, X_min_, X_max_);
    output = new TH1D(name_.c_str(), title_.c_str(), xbins.size()-1, &xbins[0]);
    return output;

  }
  TH1D* get_TH1D_from_TVectorD(string graph_title_, TVectorD *Y_values_, string Y_title_ = "", string X_title_ = "Entry #", TVectorD *Y_errors_ = nullptr) {

    auto* th1_histogram = new TH1D(graph_title_.c_str(),
                                   graph_title_.c_str(),
                                   Y_values_->GetNrows(),
                                   -0.5,
                                   Y_values_->GetNrows()-0.5
    );

    for(int i_row = 0 ; i_row < Y_values_->GetNrows() ; i_row++){
      th1_histogram->SetBinContent(i_row + 1, (*Y_values_)[i_row]);
      if(Y_errors_ != nullptr) th1_histogram->SetBinError(i_row + 1, (*Y_errors_)[i_row]);
    }

    th1_histogram->SetLineWidth(2);
    th1_histogram->SetLineColor(kBlue);
    th1_histogram->GetXaxis()->SetTitle(X_title_.c_str());
    th1_histogram->GetYaxis()->SetTitle(Y_title_.c_str());

    return th1_histogram;
  }
  TH2D* get_TH2D_log_binning(string name_, string title_, int nb_X_bins_, double X_min_, double X_max_, int nb_Y_bins_, double Y_min_, double Y_max_, string log_axis_="XY"){

    TH2D* output = nullptr;
    vector<double> xbins;
    vector<double> ybins;
    if(do_string_contains_substring(log_axis_, "X")){
      xbins = get_log_binning(nb_X_bins_, X_min_, X_max_);
    } else{
      xbins = get_linear_binning(nb_X_bins_, X_min_, X_max_);
    }
    if(do_string_contains_substring(log_axis_, "Y")){
      ybins = get_log_binning(nb_Y_bins_, Y_min_, Y_max_);
    } else{
      ybins = get_linear_binning(nb_Y_bins_, Y_min_, Y_max_);
    }

    output = new TH2D(name_.c_str(), title_.c_str(), xbins.size()-1, &xbins[0], ybins.size()-1, &ybins[0]);
    return output;

  }
  TH2D* get_TH2D_from_TMatrixD(TMatrixD *XY_values_, string graph_title_ = "", string Z_title_ = "",string Y_title_ = "Row #", string X_title_ = "Col #") {

    if(graph_title_.empty()) graph_title_ = XY_values_->GetTitle();

    auto* th2_histogram = new TH2D(graph_title_.c_str(), graph_title_.c_str(),
    XY_values_->GetNrows(), -0.5, XY_values_->GetNrows()-0.5,
    XY_values_->GetNcols(), -0.5, XY_values_->GetNcols()-0.5
  );

  for(int i_col = 0 ; i_col < XY_values_->GetNcols() ; i_col++){
    for(int j_row = 0 ; j_row < XY_values_->GetNrows() ; j_row++){
      th2_histogram->SetBinContent(i_col + 1, j_row + 1, (*XY_values_)[i_col][j_row]);
    }
  }

  th2_histogram->GetXaxis()->SetTitle(X_title_.c_str());
  th2_histogram->GetYaxis()->SetTitle(Y_title_.c_str());
  th2_histogram->GetZaxis()->SetTitle(Z_title_.c_str());

  return th2_histogram;

}
  TVectorD* get_TVectorD_from_vector(std::vector<double>& input_vector_){
    TVectorD* output = new TVectorD(input_vector_.size());
    for(int i_elm = 0 ; i_elm < int(input_vector_.size()); i_elm++){
      (*output)[i_elm] = input_vector_[i_elm];
    }
    return output;
  }
  TMatrixD* generate_correlation_matrix(TMatrixD *covariance_matrix_) {
  auto* correlation_matrix = (TMatrixD*) covariance_matrix_->Clone();
  for(int i_row = 0 ; i_row < covariance_matrix_->GetNrows() ; i_row++){
    for(int i_col = 0 ; i_col < covariance_matrix_->GetNcols() ; i_col++){
      (*correlation_matrix)[i_row][i_col] /=
      TMath::Sqrt((*covariance_matrix_)[i_row][i_row] * (*covariance_matrix_)[i_col][i_col]);
    }
  }
  return correlation_matrix;
}
  std::vector<TObject*> get_list_of_object_from_directory(TDirectory* directory_, string class_name_ = ""){
    std::vector<TObject*> output;

    for(int i_entry = 0 ; i_entry < directory_->GetListOfKeys()->GetSize() ; i_entry++){
      string object_name = directory_->GetListOfKeys()->At(i_entry)->GetName();
      TObject* obj = directory_->Get(object_name.c_str());
      if(class_name_.empty() or obj->ClassName() == class_name_){
        output.emplace_back((TObject*) obj->Clone(object_name.c_str()));
      }
    }

    return output;
  }
  std::vector<TFile*> get_list_of_opened_tfiles(){
    std::vector<TFile*> output;
    // TIter next_iter(gROOT->GetListOfGlobals());
    TList* global_obj_list = (TList*)gROOT->GetListOfGlobals();
    TGlobal *global;
    for(int i_obj=0 ; i_obj < global_obj_list->GetEntries() ; i_obj++){
      global = (TGlobal*) global_obj_list->At(i_obj);
      TString type = global->GetTypeName();
      if (type=="TFile") {
        TFile *file = (TFile*)gInterpreter->Calc(global->GetName());
        if (file && file->IsOpen()){
          // printf("%s: %s\n", global->GetName(),file->GetName());
          output.emplace_back(file);
        }
      }
    }
    // while ((global=(TGlobal*)next_iter())) {

    // }
    return output;
  }
  TMatrixDSym* convert_to_symmetric_matrix(TMatrixD* matrix_){

    auto* symmetric_matrix = (TMatrixD*) matrix_->Clone();
    auto* transposed_symmetric_matrix = new TMatrixD(*matrix_);

    transposed_symmetric_matrix->Transpose(*matrix_);
    *symmetric_matrix += *transposed_symmetric_matrix;
    for(int i_col = 0 ; i_col < matrix_->GetNcols() ; i_col++){
      for(int i_row = 0 ; i_row < matrix_->GetNrows() ; i_row++){
        (*symmetric_matrix)[i_row][i_col] /= 2.;
      }
    }

    auto* result = (TMatrixDSym*) symmetric_matrix->Clone(); // Convert to TMatrixDSym

    delete transposed_symmetric_matrix;
    delete symmetric_matrix;

    return result;

  }
  map<string, TMatrixD*> SVD_matrix_inversion(TMatrixD *matrix_){

    map<string, TMatrixD*> results_handler;

    results_handler["inverse_covariance_matrix"] = new TMatrixD(matrix_->GetNrows(), matrix_->GetNcols());
    results_handler["projector"] = new TMatrixD(matrix_->GetNrows(), matrix_->GetNcols());

    // Covariance matrices are symetric :
    auto* symmetric_matrix = TToolBox::convert_to_symmetric_matrix(matrix_);
    auto* Eigen_matrix_decomposer = new TMatrixDSymEigen(*symmetric_matrix);
    auto* Eigen_values = &(Eigen_matrix_decomposer->GetEigenValues());
    auto* Eigen_vectors = &(Eigen_matrix_decomposer->GetEigenVectors());

    for(int i_eigen_value = 0 ; i_eigen_value < matrix_->GetNcols() ; i_eigen_value++){

      if((*Eigen_values)[i_eigen_value] > 1E-5){

        for(int i_dof = 0 ; i_dof < matrix_->GetNrows() ; i_dof++){
          for(int j_dof = 0 ; j_dof < matrix_->GetNrows() ; j_dof++){
            (*results_handler["inverse_covariance_matrix"])[i_dof][j_dof] +=
                (1./(*Eigen_values)[i_eigen_value])
                *(*Eigen_vectors)[i_dof][i_eigen_value]
                *(*Eigen_vectors)[j_dof][i_eigen_value];
            (*results_handler["projector"])[i_dof][j_dof] +=
                 (*Eigen_vectors)[i_dof][i_eigen_value]
                *(*Eigen_vectors)[j_dof][i_eigen_value];
          }
        }
      } else {
        if(verbosity_level >= 4) cout << ALERT << "Skipping i_eigen_value = " << (*Eigen_values)[i_eigen_value] << endl;
      }

    }

    // No memory leak ? : CHECKED
    delete Eigen_matrix_decomposer;
    delete symmetric_matrix;

    return results_handler;

  }
  std::vector<double> get_eigen_values(TMatrixD *matrix_){
    auto* symmetric_matrix = TToolBox::convert_to_symmetric_matrix(matrix_);
    auto* Eigen_matrix_decomposer = new TMatrixDSymEigen(*symmetric_matrix);
    auto* Eigen_values = &(Eigen_matrix_decomposer->GetEigenValues());

    std::vector<double> output;
    for(int i_dim = 0 ; i_dim < matrix_->GetNcols() ; i_dim++){
      output.emplace_back((*Eigen_values)[i_dim]);
    }
    std::sort(output.begin(), output.end(), std::greater<double>());
    return output;
  }

  void save_canvas(TCanvas *canvas_, string file_name_, string sub_folder_ = "") {

    cout << WARNING << "Saving canvas in folder $FIGURES_DIR/" << sub_folder_ << endl;

    vector<string> extensions;
    extensions.emplace_back(".pdf");
    extensions.emplace_back(".png");
    extensions.emplace_back(".root");
    extensions.emplace_back(".C");

    auto old_verbosity = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;
    if(not sub_folder_.empty()) std::system(("mkdir -p ${FIGURES_DIR}/" + sub_folder_).c_str());
    for(int i_ext = 0 ; i_ext < int(extensions.size()) ; i_ext++){
      cout << WARNING << "Saving as : " << file_name_ << extensions[i_ext] << endl;
      stringstream outpath;
      outpath << "${FIGURES_DIR}/" << sub_folder_ << "/" << file_name_ << extensions[i_ext];
      canvas_->SaveAs(outpath.str().c_str());
    }
    gErrorIgnoreLevel = old_verbosity;

  }

}

namespace T2KToolBox{

  int getBeamMode(TFile* highlandTFile_){

    TTree* flattree = (TTree*) highlandTFile_->Get("flattree");

    if(flattree == 0){
      cerr << "flattree could not be found in " << highlandTFile_->GetName() << endl;
      exit(1);
    }

    flattree->GetEntry(0);
    double run = flattree->GetLeaf("sRun")->GetValue(0);

    int beam_mode = 0;
    if(int(run/1E7) == 9){
      beam_mode = 1;
      // cout << "Numu beam mode detected." << endl;
    }
    else if(int(run/1E7) == 8){
      beam_mode = -1;
      // cout << "AntiNumu beam mode detected." << endl;
    }
    else {
      beam_mode = 0;
      // cout << "Unknown beam mode : run=" << run << ". Conventions may have changed." << endl;
      // cout << "Please check : https://www.t2k.org/nd280/datacomp/production006/mcp/mcProdSummary" << endl;
    }

    return beam_mode;

  }

}
