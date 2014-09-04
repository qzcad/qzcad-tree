#include "qtscriptfunctions.h"
#include <math.h>

QScriptValue approx(QScriptContext *ctx, QScriptEngine *eng)
{
    if (ctx->argumentCount() != 2)
        return ctx->throwError("approx() takes exactly two arguments");
    if (!ctx->argument(0).isNumber())
        return ctx->throwError(QScriptContext::TypeError, "approx(): first argument is not a number");
    if (!ctx->argument(1).isNumber())
        return ctx->throwError(QScriptContext::TypeError, "approx(): second argument is not a number");
    double a = ctx->argument(0).toNumber();
    double b = ctx->argument(1).toNumber();
    const double eps = 1.0E-6;
    return fabs(a - b) < eps;
}
