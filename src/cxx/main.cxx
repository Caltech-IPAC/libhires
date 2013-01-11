#include <vector>
#include <map>

#include <log4cxx/logger.h>
#include <log4cxx/fileappender.h>
#include <log4cxx/consoleappender.h>
#include <log4cxx/patternlayout.h>

#include "Params.hxx"
#include "Exception.hxx"

log4cxx::LoggerPtr logger(log4cxx::Logger::getRootLogger());

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
      hires::Params params(argv[1],argv[2],argv[3],args);

      log4cxx::PatternLayoutPtr layout(new log4cxx::PatternLayout("\%m"));
      log4cxx::FileAppenderPtr
        file_appender(new log4cxx::FileAppender(layout,params.log_filename,false));
      log4cxx::ConsoleAppenderPtr
        console_appender(new log4cxx::ConsoleAppender(layout));
      logger->addAppender(file_appender);
      logger->addAppender(console_appender);

      logger->setLevel(log4cxx::Level::getError());

      LOG4CXX_INFO(logger,"HIRES invoked as: ");
      for(int i=0;i<argc;++i)
        LOG4CXX_INFO(logger,argv[i] << " ");
      LOG4CXX_INFO(logger,"\n");

      LOG4CXX_INFO(logger, params);

      params.compute_images();

      LOG4CXX_INFO(logger,"End Processing\n");
    }
  catch(hires::Exception &e)
    {
      std::cerr << e.what() << "\n";
      exit(1);
    }
}
