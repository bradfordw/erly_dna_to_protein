-module(dna_translate).
-behaviour(gen_server).

-include_lib("eunit/include/eunit.hrl").

-define(SGC, <<"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG">>).
-define(SGSIZE, 64).

-record(dna_ctx, {}).

-export([start_link/0, transcribe/1, codons/1, codons/2, translate/1, c2i/1]).
-export([init/1, handle_call/3, handle_cast/2, handle_info/2, terminate/2, code_change/3]).

%% API
start_link() ->
  gen_server:start_link({local, ?MODULE}, ?MODULE, [], []).

transcribe(Fragment) ->
  gen_server:call(?MODULE, {translate, Fragment}).

%% gen_server

init([]) ->
	{ok, #dna_ctx{}}.

handle_call({translate, Fragment}, _, State) ->
  Protein = codons(Fragment),
  {reply, {ok, Protein}, State};
handle_call(_Request, _From, State) ->
	{noreply, ok, State}.
	
handle_cast(_Msg, State) ->
	{noreply, State}.

handle_info(_Info, State) ->
	{noreply, State}.

terminate(_Reason, _State) ->
	ok.

code_change(_OldVsn, State, _Extra) ->
	{ok, State}.

%% Internal
codons(DnaFrag) ->
  codons(DnaFrag, []).
codons(<<>>, Protein) ->
  list_to_binary(Protein);
codons(<<A1:1/binary>>, Protein) ->
  codons(list_to_binary([A1, 0, 0]), Protein);
codons(<<A1:1/binary, A2:1/binary>>, Protein) ->
  codons(list_to_binary([A1, A2, 0]), Protein);
codons(<<A1:1/binary, A2:1/binary, A3:1/binary, Rest/binary>>, Protein) ->
  codons(Rest, Protein ++ translate([A1, A2, A3])).

translate([Base1, Base2, Base3]) ->
  case lists:member(-1, [c2i(Base1), c2i(Base2), c2i(Base3)]) of
    true -> [];
    false ->
      P = case ((c2i(Base1) * 16) + (c2i(Base2) * 4) + c2i(Base3)) of
        V when V < 0 -> binary:at(?SGC, ?SGSIZE + V);
        V -> binary:at(?SGC, V)
      end,
     [<<P>>] 
  end.

c2i(Base) ->
  case Base of
    <<$T>> -> 0;
    <<$t>> -> 0;
    <<$C>> -> 1;
    <<$c>> -> 1;
    <<$A>> -> 2;
    <<$a>> -> 2;
    <<$G>> -> 3;
    <<$g>> -> 3;
    _ -> -1
  end.

%% Tests
-ifdef(TEST).
c2i_test() ->
  0 = c2i(<<$T>>),
  0 = c2i(<<$t>>),
  1 = c2i(<<$C>>),
  1 = c2i(<<$c>>),
  2 = c2i(<<$A>>),
  2 = c2i(<<$a>>),
  3 = c2i(<<$G>>),
  3 = c2i(<<$g>>),
  -1 = c2i(<<$F>>).

translate_test() ->
  [] = translate([-1,<<$a>>,<<$t>>]),
  [] = translate([<<$t>>,-1,<<$a>>]),
  [] = translate([<<$g>>,<<$c>>,-1]),
  [<<"M">>] = translate([<<$A>>,<<$T>>,<<$G>>]).

codons_test() ->
  DNA = <<"ATGATGATAGATAGATATAGTAGATATGATCGTCAGCCATAC">>,
  Protein = <<"MMIDRYSRYDRQPY">>,
  Protein = codons(DNA),
  <<>> = codons(<<"B">>).

-endif.