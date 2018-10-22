struct task_msg_t {
	int id;
	double start_time;
	double completion_time;
	int hostname_length;
	char hostname[];
};

struct result_t {
	int size;
	char payload[];
};