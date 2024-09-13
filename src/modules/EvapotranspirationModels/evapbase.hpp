class evapT_base
{
public:
    virtual ~evapT_base() = default;

    virtual double CalcEvapT(param_base& baseparam, var_base& basevar) const = 0;

}

struct param_base
{

}

struct var_base
{

}

struct modeloutput
{
    double ET;
}
