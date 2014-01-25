
#include "CLPerf.hpp"

CLPerf::CLPerf(CLPerfMode _mode) {
    mode = _mode;
}

CLPerf::~CLPerf() {}

void CLPerf::send_msg(CLPerfCategory cat, string msg) {
    switch(cat) {
        case PERF:
            if(mode == PRINT_PERF) {
                msgs.insert(msgs.begin(), msg);
            }
            break;
        case DEBUG:
            if(mode != PRINT_ERR) {
                msgs.insert(msgs.begin(), msg);
            }
            break;
        case ERROR:
            msgs.insert(msgs.begin(), msg);
            break;
    }
}

bool CLPerf::has_msg() {
    return (msgs.size() > 0);
}

string CLPerf::recv_msg() {
    return msgs.pop_back();
}
