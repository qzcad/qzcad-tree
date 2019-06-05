#include "consoleprogress.h"


ConsoleProgress::ConsoleProgress(unsigned long expectedCount,
                                 ostream &stream,
                                 const string &s1,
                                 const string &s2,
                                 const string &s3, char spacer) :
    stream_(stream), s1_(s1), s2_(s2), s3_(s3), spacer_(spacer)
{
    restart(expectedCount);
}

void ConsoleProgress::restart(unsigned long expectedCount)
{
    count_ = nextTicCount_ = 0;
    tic_ = 0;
    expectedCount_ = expectedCount;
    stream_ << s1_ << "0%   10   20   30   40   50   60   70   80   90   100%\n"
               << s2_ << "|----|----|----|----|----|----|----|----|----|----|"
                  << endl // необходимо для вывода из буфера в поток
                     << s3_;
    if (!expectedCount_)
        expectedCount_ = 1; // упреждается деление на ноль
}

unsigned long ConsoleProgress::operator+=(unsigned long increment)
{
    inc(increment);
}

unsigned long ConsoleProgress::operator++()
{
    return inc(1);
}

unsigned long ConsoleProgress::count() const
{
    return count_;
}

unsigned long ConsoleProgress::expectedCount() const
{
    return expectedCount_;
}

bool ConsoleProgress::isExpectedCount() const
{
    return count_ == expectedCount_;
}

unsigned long ConsoleProgress::inc(unsigned long increment)
{
    if ( (count_ += increment) >= nextTicCount_ )
    {
        displayTic();
    }
    return count_;
}

void ConsoleProgress::displayTic()
{
    // use of floating point ensures that both large and small counts
    // work correctly.  static_cast<>() is also used several places
    // to suppress spurious compiler warnings.
    unsigned int tics_needed = static_cast<unsigned int>( (static_cast<double>(count_) / expectedCount_) * 50.0 );
    do
    {
        stream_ << spacer_ << flush;
    } while ( ++tic_ < tics_needed );
    nextTicCount_ = static_cast<unsigned long>((tic_ / 50.0) * expectedCount_);
    if ( count_ == expectedCount_ )
    {
        if ( tic_ < 51 ) stream_ << spacer_;
        stream_ << endl;
    }
}
