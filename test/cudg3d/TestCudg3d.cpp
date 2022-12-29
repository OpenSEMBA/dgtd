#include "gtest/gtest.h"

#include "cudg3d/Cudg3d.h"
#include "parsers/json/Parser.h"

using namespace SEMBA;
using namespace SEMBA::dgtd;

class SolverTest : public ::testing::Test {
public:
    FileSystem::Project getProjectFilename(const std::string& project)
    {
        FileSystem::Project res{ "./testData/" + project + "/" + project + ".smb.json" };
        EXPECT_TRUE(res.canOpen()) << res;
        return res;
    }

    UnstructuredProblemDescription readProject(const std::string& name)
    {
        return Parsers::JSON::Parser(getProjectFilename(name)).read();
    }
};

TEST_F(SolverTest, resonant_cube) 
{ 
    auto smb{ readProject("resonant_cube") };
    Cudg3d::Options opts;
    ASSERT_NO_THROW(Cudg3d cudg3d( smb, opts ));
}