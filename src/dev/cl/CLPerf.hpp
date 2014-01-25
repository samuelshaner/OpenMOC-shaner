#include <string>
#include <vector>

using namespace std;

enum CLPerfMode {
    PRINT_PERF,
    PRINT_DEBUG,
    PRINT_ERR
}

enum CLPerfCategory {
    PERF,
    DEBUG,
    ERROR
}

class CLPerf {
private:
    vector<string> msgs;
    CLPerfMode mode;

public:
    CLPerf(CLPerfMode);
    ~CLPerf();

    void send_msg(CLPerfCategory, string);
    bool has_msg();
    string recv_msg();
}
