void systemVerbose(const string call)
{
//     if (call[0] == 'r' && call[1] == 'm')
//     {
//         Info << "Not executing " << call << endl;
//     }
//     else
//     {
//         Info << call << endl;
        bool err = system(call.c_str());
        if (err)
        {
            FatalErrorIn(call) << exit(FatalError);
        }
//     }
}

static string systemCall;
static string systemCall2;

