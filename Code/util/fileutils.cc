#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <unistd.h>
#include <dirent.h>

#include <sys/dir.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "log/Logger.h"
#include "util/fileutils.h"

namespace hemelb
{
  namespace util
  {
    namespace
    {
      // Returns true if the file with the given name exists for reading,
      // false otherwise.
      //bool file_exists(const char * filename);

      // Function to select directory contents that are not "." or ".."
      // int selectOnlyContents (direct_t *entry);


      // Return true if file exists for reading, false if not.
      bool file_exists(const char * filename)
      {

        if (access(filename, R_OK) == -1)
        {
          return false;
        }
        return true;
      }

    }

    // Check the existence of a critical file - exit if it's not there
    void check_file(const char * filename)
    {
      if (!file_exists(filename))
      {
        log::Logger::Log<log::Info, log::OnePerCore>("Cannot open file %s\nExiting.", filename);
        exit(0);
      }
    }

    // Function to select directory contents that are not "." or ".."
    int selectOnlyContents(direct_t *entry)
    {
      if ( (strcmp(entry->d_name, ".") == 0) || (strcmp(entry->d_name, "..") == 0))
      {
        return 0;

      }
      else
      {
        return 1;
      }

    }

    void ChangeDirectory(const char * target){
      chdir(target);
    }

    void ChangeDirectory(const std::string & target){
      chdir(target.c_str());
    }

    void GetCurrentDir(char * result, int bufflength){
      getcwd(result,bufflength);
    }

    std::string GetCurrentDir(){
      char buff[1000];
      GetCurrentDir(buff,1000);
      return std::string(buff); // return by copy.
    }

    // This copied from BOOST. TODO: Use boost
    const char * GetTemporaryDir(){
      const char *dirname;
      dirname = std::getenv("TMP");
      if(NULL == dirname)
        dirname = std::getenv("TMPDIR");
      if(NULL == dirname)
        dirname = std::getenv("TEMP");
      if(NULL == dirname){
        //assert(false); // no temp directory found
        dirname = ".";
      }
      return dirname;
    }

    // Delete all files within a directory.
    int DeleteDirContents(std::string pathname)
    {
      struct direct **files;

      int file_count = scandir(pathname.c_str(), &files, selectOnlyContents, alphasort);

      char filename[1024];

      for (int i = 0; i < file_count; i++)
      {
        snprintf(filename, 1024, "%s/%s", pathname.c_str(), files[i]->d_name);
        unlink(filename);
      }
      return 0;
    }

    // Detect whether a directory exists.
    bool DoesDirectoryExist(const char *pathname)
    {
      struct stat st;
      return stat(pathname, &st) == 0;
    }

    bool FileCopy(const char* iOriginalPath, const char* iNewPath)
    {
      std::ifstream lSource;
      std::ofstream lDestination;

      // open in binary to prevent jargon at the end of the buffer
      lSource.open(iOriginalPath, std::ios::binary);
      lDestination.open(iNewPath, std::ios::binary);

      if (!lSource.is_open() || !lDestination.is_open())
      {
        return false;
      }

      lDestination << lSource.rdbuf();

      lDestination.close();
      lSource.close();

      return true;
    }

    // Function to create the directory of given path, which user group and anyone
    // can read write and execute.
    void MakeDirAllRXW(std::string &dirPath)
    {
      mkdir(dirPath.c_str(), 0777);
    }

    std::string NormalizePathRelativeToPath(std::string inPath, std::string basePath)
    {
      // If it's an absolute path, just return it
      if (inPath[0] == '/')
      {
        return inPath;
      }

      // Going to check if it's a directory
      std::string baseDir;
      struct stat st;
      stat(basePath.c_str(), &st);
      // Assume it's a regular file in case it doesn't exist
      st.st_mode = S_IFREG;

      if (st.st_mode == S_IFDIR)
      {
        // It's a directory
        baseDir = basePath;
      }
      else
      {
        // Not a dir, find the last slash
        unsigned long lastSlash = basePath.rfind('/');
        if (lastSlash == basePath.npos)
        {
          // No slashes, so the baseDir is just the working dir
          baseDir = ".";
        }
        else
        {
          // Has slashes, return up to the last
          baseDir = basePath.substr(0, lastSlash);
        }
      }

      // Make sure it ends in a slash
      if (baseDir[baseDir.size() - 1] != '/')
      {
        baseDir += "/";
      }

      //Append the path of interest
      return baseDir + inPath;
    }

  }
}
