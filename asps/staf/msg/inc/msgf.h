/*
	This contains prototypes for the msg library, standard FORTRAN calls.

	Created  2-Dec-1994   Robert Hackenburg
*/

#define TRUE -1
#define FALSE 0

	void message_( const char *msg, int *one, int *ID, int len );
	void message_out_( const char *msg, int *one, int len );
	void msg_define_class_( const char *Class, const char *State, int *CountLimit, int *AbortLimit, int Clen, int Slen );
	void msg_count_( const char *Prefix, int len );
	void msg_disable_( const char *Prefix, int len );
	void msg_display_( const char *msg, int *one, int *ID, int len );
	void msg_display_out_( const char *msg, int *one, int len );
	void msg_display_and_echo_( const char *msg, int *one, int *LUN, int *ID, int len );
	void msg_display_and_echo_out_( const char *msg, int *one, int *LUN, int len );
	void msg_enable_( const char *Prefix, int len );
	int  msg_enabled_( const char *Prefix, int *ID, int len );
	void msg_get_lun_( int *Terminal_LUN, int *Journal_LUN );
	void msg_ini_( int *Journal_LUN );
	int  msg_journal_close_( void );
	void msg_journal_off_( void );
	void msg_journal_on_( void );
	int  msg_journal_open_( const char *FileName, int len );
	void msg_journal_page_( void );
	void msg_lun_page_( int *LUN );
	void msg_mark_( const char *Prefix, int *ID, int len );
	void msg_node_name_( const char *NodeName, int len );
	void msg_nocount_( const char *Prefix, int len );
	void msg_set_abort_limit_( const char *Prefix, int *Limit, int len );
	void msg_set_by_command_( const char *Command, int len );
	void msg_set_from_file( int *LUN );
	void msg_set_limit_( const char *Prefix, int *Limit, int len );
	void msg_set_lun_( int *TERMINAL_LUN, int *JOURNAL_LUN );
	void msg_set_summary_mode_aborted_( int *Mode );
	void msg_set_summary_mode_active_( int *Mode );
	void msg_set_summary_mode_counting_( int *Mode );
	void msg_set_summary_mode_inactive_( int *Mode );
	void msg_set_summary_page_length_( int *Page_Length );
	void msg_set_timestamp_cpu_( int *Mode );
	void msg_sort_( void );
	void msg_summary_( int *LUN );
	void msg_summary_cpu_( int *LUN );
	void msg_summary_event_( int *LUN, int *EVENTS );
	void msg_time_stamp_( int *LUN );
	void msg_time_stamp_out_( int *LUN );
	void msg_to_journal_( const char *msg, int *one, int *ID, int len );
	void msg_to_journal_out_( const char *msg, int *one, int len );
	void msg_to_lun_( const char *msg, int *one, int *LUN, int *ID, int len );
	void msg_to_lun_out_( const char *msg, int *one, int *LUN, int len );
