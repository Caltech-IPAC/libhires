#include <vector>
#include <map>

#include "Hires.hxx"
#include "Exception.hxx"

int main(int argc, char* argv[])
{
  if(argc<5)
    {
      std::cerr << "Usage: hires data_type IN_prefix OUT_prefix"
        " param_file1 param_file2 ...\n";
      exit(1);
    }

  std::vector<std::string> args;
  for(int i=4;i<argc;++i)
    args.emplace_back(argv[i]);

  try
    {
      hires::Hires params(argv[1],argv[2],argv[3],args);
      params.compute_images();
    }
  catch(hires::Exception &e)
    {
      std::cerr << e.what() << "\n";
      exit(1);
    }
}
