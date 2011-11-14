#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/XmlOutputter.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>
#include <cppunit/TextTestProgressListener.h>
#include <stdexcept>
#include "unittests/lbtests/lbtests.h"
#include "unittests/vistests/vistests.h"

int main(int argc, char **argv)
{
  std::string testPath = (argc > 1)
    ? std::string(argv[1])
    : "";
  // Create the event manager and test controller
  CppUnit::TestResult controller;

  // Add a listener that colllects test result
  CppUnit::TestResultCollector result;
  controller.addListener(&result);

  // Add a listener that print dots to stdout as test run.
  CppUnit::TextTestProgressListener progress;
  controller.addListener(&progress);

  // Add the top suite to the test runner
  CppUnit::TestRunner runner;

  runner.addTest(new hemelb::unittests::lbtests::LbTestSuite());
  runner.addTest(new hemelb::unittests::vistests::VisTestSuite());

  try
  {
    std::cout << "Running " << testPath;
    runner.run(controller, testPath);

    // Print test XML output to stderr
    CppUnit::XmlOutputter outputter(&result, std::cerr);
    outputter.write();
  }
  catch (std::invalid_argument &e) // Test path not resolved
  {
    std::cerr << std::endl << "ERROR: " << e.what() << std::endl;
    return 1;
  }

  return result.wasSuccessful()
    ? 0
    : 1;
}
