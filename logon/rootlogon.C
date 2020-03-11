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

}

#include <sys/types.h>
#include <sys/stat.h>

namespace TToolBox {

  // Forward declarations
  std::vector<std::string> split_string(std::string input_string_, std::string delimiter_);
  std::string join_vector_string(std::vector<std::string> string_list_, std::string delimiter_, int begin_index_, int end_index_);

  ////////////////////////////
  // Src
  static double last_displayed_value = -1;
  static int verbosity_level = 0;
  static Color_t colorCells[10] = {kOrange+1, kGreen-3, kTeal+3, kAzure+7, kCyan-2, kBlue-7, kBlue+2, kOrange+9, kRed+2, kPink+9};

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
  std::string join_vector_string(std::vector<std::string> string_list_, std::string delimiter_, int begin_index_, int end_index_) {

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


}
