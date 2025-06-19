template<class data>
class base_step
{
public:
	base_step(data& _d) : d(_d) {};
	virtual ~base_step() = default;

	virtual void execute() = 0;
protected:
	data& d;
};
