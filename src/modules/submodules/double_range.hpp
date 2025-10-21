namespace double_range
{
    class Unit // Real number on the interval [0,1]
    {
        double _c;
    public:
        Unit(const double c);
        ~Unit() = default;
        
        operator double() const;
    };
};
