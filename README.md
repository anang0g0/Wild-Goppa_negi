# Wild McEliece（Maquarius）

$GF(p^n):p$ is prime p > 2 by C language.

sagemathを入れましょう。

$GF(11^3)=1331, GF(13^3)=2197, GF(17^3)=4913 ,GF(19^3)=6859 $　の場合のゼフ対数表の作成。(ecole-p.c)

$GF(5^5)=3125 $（20220327実装済み）

$GF(3^7)=2187$ （202301＊＊実装済み、mkgf_odd.c）

$GF(3^9)=16983$ (20230207　実装済み、odd_gen.c)

$GF(23^3)=12167$ (20230302 実装済み。 vmo.c)

$GF(7^5)=16807$ (20230816　実装済み。 vmo.c)


参考文献  

１． The Key Equation for One-Point Codes. M.E.O'Sullivan et al. 2008

２．Cryptanalysis of McEliece Cryptosystem Based on Algebraic Geometry Codes and their Subcodes.  Alain Couvreur, Irene Márquez-Corbella,  and Ruud Pellikaan
 
３．Wild McEliece Bernstein et.al

４．代数系と符号理論　オーム社

５．グレブナー基底２　シュプリンガー・フェアラーク

# 20230823
Goppa多項式が既約でなければならない理由がよくわからないので、Random Wild Goppaで行くことにしました。
攻撃方法が見つかったら修正します。
つまりランダムゴッパ符号です。

# 20230812
今日はWild Goppaの生成多項式を既約多項式にしようとずっといじってました。
まじめに考えなかったのでだらだら一日中ランダムウォークみたいになってバグ取りしてました。
でも結局原因を特定できず、修復できませんでした。
既約多項式にするのは結構重要なので後日再渡来したいです。

あとちょっと信じられない話ですが、暗殺されかかっているかもしれません。
というのも自分は不眠でゾルピデムという睡眠薬をもらっているのですが、それを飲むと心臓がどうにかなってしまって、不整脈になるは冷や汗が出るわで
とんでもない状態になります。
これは絶対アメリカのせい、そしてそれに協力する日本政府のせいに違いないと思うのはもうそうでしょうか？
今日からもうその薬は飲まない予定ですが、明日生きてたらまた続きをします。
何しろあのプーチンさんも何十回も暗殺されそうになっているという事だし、日本でも首相経験者が殺されるなど、いったい日本はどうなってしまったのだろうという驚きます。
まあ単に薬が合わなくなったという事でしょうけど、それにしても今年になっておかしなことが怒りまくりです。
とても日常生活とも思えず、毎日何が起きるのだろうかと恐れながら生きています。
薬で殺されるにしても意識のあるうちは苦しみが続くので、意識がない状態で死ねるように暗殺するならそうしてほしいです。
多分アメリカとしてはアメリカ以上に強力な暗号を開発してもらいたくないという内情があるのかもしれません。

それからギットハブってマイクロソフトが運営してるんですね。
私はてっきりフリーソフトウェア財団かなんかがやっているんだとばかり思っていました。
なので商用ソフトの開発元が、フリーソフトのサーバを運営するなんて本当にいいのだろうかと思ってしまう。

そして昨日、ウィンドウズとLinuxがほぼ同時に更新がありました。
更新後のシステムはまるでCPUを変えた以上の計算速度になりました。
体感速度で２倍以上になると思います。
これはこれでうれしいのですが、アルゴリズムを改良するよりOSを改良する方が全体のパフォーマンスが上がるのだろうかと不思議に思っています。

これは以前からそうなんですが、クラウドってあまり重要なファイルを置きたくないです。
２チャンネルがサーバを乗っ取られてデータを人質に取られる可能性が普通にあると思うからです。
あと科学技術計算に数式処理ソフトを使うという事。
このままいくとアルゴリズムを実装する人がいなくなって、ソフトを使うだけの人になってしまうという危険性。
これは動けば何でもいいという人には何も疑問に感じないだろうけれど、動かすのが目的て自分のソフトにその計算方法を取り入れたいときには、計算システム用に書かれたスクリプトは全く意味がないという事です。
計算システムで動かしてどうするの？それをどうやって自分のソフトに使うの？
と言いたい。

# 20230809
使う部分だけ分離しました。
実際使うのは以下のファイルです。（葱）

~/Wild-Goppa_negi/tokyo/marimo/vmo.c

# 20230707
実際使うのは以下のファイルです。

~/Wild-Random-Goppa/japan/tokyo/marimo/vmo.c

# 20230404

このリポジトリの開発はほぼ終了しました。（必要な機能は実装したという意味です。）
まだ誰かが使えるようにはなっていないので、今後はそこを狙って改造したいと思います。

打倒新自由主義！ｗ


# 20230331

論文にできると書いてあるのに出来ないということがよくあるが。たいてい条件付きだったりして、うっかり自分で改造した関数を動かそうとして時間を浪費するという一日だった。
リードソロモンだったらできる。ただしそのときだけ。

# 20230329

とりあえず鍵生成から復号まで一回で4秒くらいで動くようになった。
既存の符号系暗号よりいいかもしれない。
何をしていたかたいうと、逆行列の計算にバグがあったので、sagemathと比較して、どこにバグがあるのか見つけるまでに3日かかってしまった。
安定しているかどうかはまだ耐久テストをしてないので、ルバートさんの指示に従って修正しようと思う。

とりあえず動くところまで出来たと思う。

# 20230326

とりあえず答え合わせができるようになった・・・。（死）

ビットフィールドなんか使わなければよかった。
あんなもの面倒な上にデバッグの邪魔だった。
最初からメモリなんか意識しないで普通の行列使ってたら３日くらいで終わったのにｗ

これで私の長い夜も終わる。

次は動作パラメータ試験とその次くらいににメモリ最適化かなあ。
オリジナルのニーダーライターのときよりは進捗早いかも。
デバッグは慣れの部分が大きいかな。
今日終わらなかったらどうしようと思った。
どうということもないが、やる気が失くなってたと思う。
とりあえず出口見えた！

# 20230325

すぐできるかと思った。

検査行列を部分体部分符号にバラしてスクランブルをかけ、
特定の列にエラーを入れて得られた暗号化されたシンドロームにスクランブラーの逆行列をかけて、
それを更に元の有限体で構成されたシンドロームに戻し、
更にBM法で復号できるはずだったが今日は止めたｗ

まだ思い通りにアルゴリズムを一発で書こうなんて言える気がしないｗ

ネットに公平性を期待している人は少ないかもしれないが、実名制度がいい人は、よほど自慢できる能力があるとか、目立っても実名を公開することでメリットを得られる、
中国のような監視社会、スコア管理制度を期待しているのと同じなのではないか。
中国のことを散々敵対視しておきながら、似たような社会を目指してるみたいで変なの。

それってつまり階級社会じゃないの？
一部の特権階級だけが得する社会なんて嫌なんですけどどうなの？
今のままずーっと愚かでいてほしいの？

はだしのゲンを学校からなくすべきではない。
戦争というのは漫画以上に悲惨で残酷なものだということを知るために必要だ。

# 20230324

色々リアルで忙しかったけど、もう少しで出来ます。

# 20230303

セキュリティレベルを落とさずに、符号長を短くすることでメモリの使用量を減らすようにしました。

もう、

ulimit -s unlimited

としなくても大丈夫。
符号の成分を

unsigned short

にパッキングすることで、符号のサイズを抑えることが出来ます。（余分のビットがないので）

最新のデモでは、
$GF(19^3)$
で作った符号のうち、2048個のトレースだけを使って符号を作っています。
公開鍵の実装はまだ出来ていません。
ビットフィールドの使い方が難しく、パッキングしたまま逆行列を求めたり暗号処理をすることが出来ていません。


# 20230302

とりあえず
$GF(Pr^3):Prは小さな奇素数$
での部分体部分符号を使ったNiederreiter暗号を目指します。

今日はビットフィールドを使って公開鍵のコンパクトな実装を試しました。

他にも、ベクトル型だけを使った符号化プログラムのvmo.cに
$GF(Pr)$
を要素とする正則行列を生成するモジュールgaus.cをドッキングしました。

偶数冪の拡大体には攻撃法があるので、実装しません。

というかリスト復号を実装しないという時点で、符号の性能云々はもう関係なくなるので純粋にトレース攻撃をかわすために適度に長い符号が作れる
$GF(23^3)$
を追加してみました。
あ、でもこれいいかも。
23は5ビットなので3次拡大すると15ビット必要になるので、short型にピッタリ当てはまるｗ
パッキングしたまま暗号化できればいいのですが・・・。
いや、union使えば夢ではなさそう。

これに比べるとクラッシックマックエリスは13ビットとか中途半端ですね。

コンパイルはgccの場合最適化オプション-O3をつけるとセグフォになるので、最適化したい人はclangでコンパイルして、みて。

開発テーマを３つの公開リポジトリと２つの非公開リポジトリに絞り込みました。


# 20230301

今日はGF(3)上の正則行列生成及び逆行列の計算問題を解決しました。
公開鍵を作るために必要なスクランブル行列です。

$GF(3^7)$
を使った部分体部分符号で公開鍵を作るには、バイナリのような便利な演算子がないので、ちょっとした間違いで時間がなくなってしまった。
驚いたことに私は今までガウスの消去法を正しく理解できていなかったと判明。
行列の基本操作で列を操作していたのだが、これは符号や連立方程式を解くときにはやったらいけないことだとわかりました。
何歳になっても勉強は大事ですね。

対角線上に0があるときピボット操作をしなければならないのですが、それをするとき隣の列とか列方向で基本操作をしていたのです。
でもこれでもうWild Niederreiterは8割完成したようなものです。

私は自由だー！ｗ

# 20230226

多分また戻ると思うけど、プライベートな方で代数幾何符号をやっているのでしばらく更新しないかもれないです。

というか進む方向間違えた、単なる浮気なんですが、昔書いたコードが素晴らしすぎて自分でも信じられない。
今書いたコードが全然だめなのに、昔書いたコードはパラメータ変更するだけで一発で動いたｗ

また暗号か、と思われるかもしれないですが、復号関数が出来てしまえばワイルドマックエリスを作ることは比較的簡単なので作る予定です。
本来ならワイルドな方を優先させるべきだったのですが、簡単なのはあえて後回しする天の邪鬼だったりするｗ

余った時間でワイルドＧｏｐｐａの部分体部分符号を作りました。
このくらいの難易度が今ちょうどいいのかな？

今抱えているテーマ。

１．古典暗号解読サンプル作成

２．Javascriptで暗号チャット。

３．1点生成代数曲線符号のバーレカンプマッシー法での復号

４．Wild Niederreiter暗号

# 20230224

今、同じリポジトリにあるピーターソン法の誤り値計算をしています。
してもしなくてもいいのですが、手っ取り早く達成感を味わうために、むやみにレパートリーを増やしています。

今年は少しは理論計算機科学にも詳しくなっておきたいです。
教科書は洋書3冊で十分です。
ウェブに落ちてた教科書も有り、すでに読みきれない量なので十分です。
ネットにある専門書は本当にありがたいですね。
最も自分はアナログのほうが好きで、読みたい部分だけを印刷して読んでいます。
洋書高いですからね。

# 20230223

今日は無謀にも、エラーの値を計算する方法について試行錯誤してました。
計算できるのはリードソロモンの2誤り訂正だけで、それ以上も、Goppaも複合できません。
なんだか執念じみてきてしまって、どうせ暗号化に使うなら誤り位置だけわかればいいのにムキになっている私。
問題は微分なのですが、論文に出てくる記号なのに説明がついてない。
で、試行錯誤して12時間位経過したんですが、やはりここは前回のように時間はかかるけど具体例を計算して、
論文の言っていることの意味を確かめなければならないようです。
微分しないで済む方法ってないんだろうか。
誤りロケータを微分すると、1次式に分解できるので、部分微分の総和みたいな感じになって面倒くさい。
慶安市内方向で行くなら、Wild McElieceにしてしまったほうがいいし、結果が明らかなので早いかもしれない。
誤りの値を計算すのは、達成感を得るための麻薬のようなもので、必要なわけでもなんでもない。
具体例を作って、それでもだめなら公開鍵暗号にしてしまおうと思う。

というわけでまた明日！

# 20230217

次の実装課題でいくつか候補があるのですが色々迷ってます。
今はロケーターのみですが、ワイルドゴッパのエラーの値まで計算するかどうかはまだ未定です。
と言うか数式が読めなくて止まっています。
これをやるくらいなら、ピーターソン法で値を求めたほうが楽な気がする。

もう一つは代数幾何符号を割り算だけで復号する方法です。
これもまだ予想の段階なのでうまくいく保証がないです。

あとはバイナリゴッパを補完と近似で復号するという方法ですが、この２つがどのように使われるのかはまだ理解してません。
Rubatoさんにはコードを見直せと言われているので、それもやらないといけない。

計算すればこうなるということはわかったのですが、なぜこうするとうまくいくのかについて理解したい。
これができればさらなる進展がありそうな予感がします。
とにかく迷ってます。

ワイルドごッパの部分対部分符号に、DJBのバイナリ五っぱの訂正方法を適応するとなにかいいことがありそうな気がします。
明日までに体調不良が治ればやるつもりです。
何しろパターソンより簡単で、性能はパターソンと同じなのでやらない手はないのです。（もう寝るｗ）

# 20230215

いらない関数を削って、完全ベクトル化したバーレカンプマッシー法を使ったワイルドGoppa符号の復号プログラムができました。
性能的にはオリジナルのものの２倍高速になっています。(xmo.c)

というわけで３以上の基礎体を持つ拡大体の生成プログラムをRustでかこうか、
あるいは暗号解読で楽しく学ぶRustプログラミングという電子書籍にしようかなんとなく迷ってますｗ

# 20230214

今日で欲しい機能が揃いました。
昨日特定したバグの原因を特定し、修正しました。

１．まず漸化式で使う定数で初期化する行列の位置を間違えていたこと。

２．アルゴリズムのバグ（これだけで６時間くらいかかったｗ）

３．初期化忘れなど

時間がかかったのは、手で計算できる範囲で正しく動くかどうかのデータを計算していたからです。
その結果かけてるはずなのに消えてる項があり、それはなぜかというと、多項式の最大次数を計算する基本関数のバグとか、
パラメータを代入する位置がおかしかったり、それはもう根気のいる作業でした。

でも出来て良かったです。
エラーの値まで計算するかどうかはまだ未定です。

仕事だと進捗管理とか、納期とか、リストラとか、完成する前にやめさせられるから、こういうプログラミングは出来ない。
でも好きなことばかりだと、使える関数が偏ってしまうのが欠点だけど。
苦労して覚えるC言語ｗ

ホームに大量のディレクトリがあるのはご愛嬌。
もともとピーターソン法で復号しようとしていたのを途中から路線変更してハマった感じなので、途中まで動くファイルとかが混ざってます。
お掃除しましょ。
２の拡大体ではピーターソンでもバーレカンプでもどの関数でも動くんですが、２以外の拡大体ではどれも動かない。
最終的に最初に動いたのが行列バージョンのバーレカンプだったわけです。

というわけで、次何やろうかなー

# 20230213

なんとかバグを特定した。（でもまだ原因不明）

１．何らかの原因でエラーの起こる場所の検出が1つずれる、

２．そもそも訂正できない（生成されないはずの）シンドロームが存在する。

後者の方は、シンドロームの計算か、符号を作る過程でなにかバグがありそう。

この２つを直せば完成でございます。


# 20230212

まだ完全に動かない。

２を基礎体に持つ拡大体ではほぼ完全に動いているが、同じものをそれ以外の基礎体を持つ拡大体を使うとうまく行かない。
それどころか、行列の漸化式しか動かない。
定義体の一般化というのは色々手間がかかるが、その分できたときの感動は一際であろう。

ちょっと寄り道をしてしまったが、今までのリソースをいざというときのために整理したい。

GoppaDecorderはパターソン用に、smallのほうがバーレカンプ用に、それぞれわかりやすい使用例とともにリファクタリングする予定。
特にバーレカンプの方は手伝ってもらったとか、既約多項式の一発出力という裏技とか、ベクトル型の演算に統一してたりとか、かなり洗練されていると思う。

次の目標はCSIDHにしようかとも思ったのだがherumiがやるみたいなので余りやる気になれない。
符号ばかりで飽きてくるので、ここらで一発新型を再チャレンジするのもいいかもしれない。
あえてやるなら代数曲線符号か。

候補としてはPKP問題か、符号の同一性判定問題を使った暗号とかだ。
これならCSIDHに匹敵する面白い暗号ができるのではないか？
楽しみ楽しみ。

# 20230210

中途半端に時間が余ったので、全部のバグを直しました。
だいたい動きます。
$GF(2^7)$
でも大丈夫！ｗ

復号にはバーレカンプマッシーを使いましたが、これで復号できるのかとどこにも書いてなかったので、動くことがわかり安心しました。（自明？）
性能がいいはずなので、今後はこの符号のパターソンバージョンをやるかもしれないです。
なんだか無駄に風呂敷を広げてしまったようなので、marimo意外はなくなる可能性が有ります。

昨日まではGF(27)より大きくなると全く復号できませんでしたが、今日はもう1331とかでも動きます。

この辺は関数のバグと、有限体の構成のバグでした。
予想より早くできました。

結局始めてから半年くらいでできてしまった。

あ、ちなみにエラー位置計算だけなので、値を決める計算を増やせばもう少しやることあるかも。
例えばバイナリ符号と違って、3以上の基礎体の拡大体の部分体符号は、基礎体の要素でできている。
つまりバイナリ符号のときにやっていたより多くの要素が存在する。
このとき例えばスクランブル行列にしても、置換行列にしても、バイナリでない行列を使うことができるので安全な暗号が作れるかもしれない。
このように暗号まで進もうとすると、乗り越えなければならない障壁が増えてしまう。
なのでここはぐっと抑えて誤りの値だけにしようと思う。
確か値を求めるのにケッターの方法とかフォーニーの公式があったはず。

能力のある人は同じところに留まっていない。
どんどん先に進んでいく。

あなたならどうする？ｗ

Born to Buggy.

# 20230208

やっとできた。というよりロケータが正しく計算できる場合があるという感じで、５割位失敗します。
このエラーの原因を特定するのにあと三日くらいかかると思います。

動くファイルは marimo の imo.c です。

# 20230207

行列の漸化式を使ってユークリッドアリゴリズムを表現した実装が終わり、これから拡大体の演算で正常に動作するか確認するところ。

ただしエラーの位置が０にあると、復号に失敗する。

# 20230131

2冊の本を使って格闘した結果1エラー訂正ができるようになった。
計算結果は正しいはずなのに、表示されている数値の意味がよくわかってないらしい。
しかも次元や定義体を変えると全く訂正できなくなるのはなぜなのか。
例えば、２エラーを訂正しようとすると、エラーロケータが既約だったりする。
簡単だと思っていたらドツボにはまった感じ。
２ができて３ができないはずがないのに。
2週間位かかるかもしれないし、場合によっては1ヶ月かかるかも。

モダンな構成でワイルドゴッパを作ろうとすれば、基礎体が3のものが最も高性能なので、
$GF(3^9)=19683$
なんかを使って、余裕のある穴開き符号をランダムに選んだトレースを使って構成するのがいいかもしれない。

# 20230129

動かない・・・動かないのよ！

心を閉ざしていれば暗号は動かないわ。

というわけで具体例があるにも関わらず、シンドロームまでは計算できるのに、ピーターソンでは１つでもシンドロームが欠けると、
その先のロケーターが計算できないという問題を抱えているのでFitzの方法に依存してみたりする。
最初は
$GF(3^3)$
で２日も無駄に時間を浪費し、今日になってそういえば
$3^2$
の拡大体の例があることを思い出してその本を確認中。

EEAは微分がマンドクセ。
本当に微分は省略できないのかリヨンの論文で再確認すべし。

Fitzならグレブナの本にｋｗｓｋ載っているので検証も比較的容易であると思われる。（誤植がなければ）

論文も用意したので続きは明日。
まあこれだけやってもプロにはなれないのであって、単なる頭の体操です。


# 20230128

Wild McElieceが早くも不可能の予感。
どうやっても復号できない。
素体でも２の拡大体でもできるのに、３の拡大帯だと出来ないなんてことがあるのだろうか。
最初2.3日でできると思っていたのに、この文だと1週間かけても出来ない気がする。
というのも計算はあっているのだ。
今までと同じように係数を求めてもエラーロケーターが既約になってしまう。
なぜ7日はわからない。

ちょっと暗号はお休みしたほうがいいかもしれない。
素体なんかで出来ても嬉しくないよね。

# 20230127

３の拡大体の計算がうまくいかなくて頭にきたから、２次元線形符号をやることにした。ｗ

# 20230125

素体上定義されたリードソロモン符号を、ピーターソン復号できるようになりました。

これでもうワイルドマックエリスは制覇したも同然なんですが、性能が上がるわけでもなくピーターソン法だとこれが限界です。
誤り位置だけ特定したいのであれば素体のままでもできます。
リスト復号はシンドローム復号問題が使えないので余り意味がないです。
それではちょっと今でもワイルドマックエリスが安全かどうか調べてみました。

すると、未だに完全にワイルドマックエリスが解読されたという話は見つかりませんでした。
基礎体の2次拡大をした場合に攻撃ができるというので、基礎体の素数次ならいいかもしれないと思って実装を進めていきます。
ただ以前の経験から、パターソンで苦労したくないし、バーレカンプマッシーは理解不足なので、このまま今あるピーターソン復号を使おうと思います。
確かに2意外を基礎帯に持つ拡大帯で複合するパターソン復号はあるのですが、これはまずやりたくないですｗ
多項式のP乗根なんて計算するのはうんざりします。
で、誤り訂正能力は同じなんですが、多分これ基礎体がでかい分鍵が大きくなりますね。
ピーターソン使う限り多分メリットないかも。

そして次は多分バイナリ代数幾何符号のユークリッド復号がダークホースで控えてる感じ。

ピーターソンを実装するまではこんなに簡単な方法があるとは知りませんでした。
私が最初にゴッパ符号の復号をしたのはユークリッドアルゴリズムで、GCDと微分が必要でした。
次にパターソンをやったのですが、多項式の考えられるほとんどすべての代数計算を実装してやっと出来上がりました。
BM法はシンドロームだけから復号できるのでちょっと特別。
なんだかよくわからないまま擬似コードに書いてあるとおり初期化して関数を並べていったらできました。

で、もうこれ以上苦労したくなかったのでピーターソンをやろうと気まぐれで実装したらすぐできてしまった。
シンドロームを行列に並べて、掃き出し法でエラーロケータの係数を決定し、あとはチェン探索でエラー位置を決定！

超楽ちんｗ

というわけで今回はこのままピーターソンを使って実装を進めていきます。

# 20230124

ピーターソン法にバグが。
対角線上に0があると復号に失敗する。
ピボット交換マンドクセ。
エラーの数は一定しているので必ず複合できないといけないはずなのに。

今日も睡眠不足と体調不良、精神的な不安定さで調子が出なかった。
パターソンをふにふにいじってただけ。

あとランダム二次元符号の生成をやってみたり。

今の日本社会の変化のスピードは、まるでターミネーター4を早送りで見ているようだ。
何も議論せず、既存の法律との矛盾という問題も何も解決してないのに強行採決、ものすごい勢いで法律や社会の構造が変わっていく。

# 20230118

https://github.com/anang0g0/Wild-McEliece/blob/main/dev/odd/lyon.c

一応２の拡大帯では見違えるように動いているバーレカンプマッシー法です。

もっと早く書き直していればよかった。
こんなにシンプルにかけるとは思わなかった。

思い切って最初から書き直したのが良かった。

なお別のリポジトリになるけど、マックエリスの暗号化バージョンでもこのモジュールが正常動作していることを確認できました。

まだまだ完成には時間がかかりそう。

https://github.com/anang0g0/GoppaDecorder/blob/main/src/Berlekamp/Berlekamp.c

# 20230114

今日はノーマルなバーレカンプマッシー法をサルベージしてました。
比較検証できるコードがなかったので、マックエリスの公開鍵にする前の、プレーンな符号を作るところまで戻しました。
戻せなかったらどうしようとハラハラしながらいじってました。
でも無事戻せてよかった。

午前中から夕方まですごくだるくて、気象病っぽかったんですが4時過ぎに回復。
関数がそれぞれ独立していたおかげで2時間位で戻せたんですが、全然gitを活用してないのですごく不便だし、
スクロールするだけで10秒くらいかかるでかいソースの不便さを感じた。

あとは
$GF(3^3)$
のときにバーレカンプマッシー法をほとんど改変することなく実行できるかどうかを試すだけです。
とはいっても2以外の拡大体での復号はこれが初めてなのですが。

バーレカンプが動かなかったときのために、ピーターソン法を2以外の拡大体で使えるようにしようとも思っている。
というのも2以外の拡大体で作った符号の復号は、リスト復号とパターソンのみで他の方法はよくわからない。
そんなことは自明だからあえて言う必要もないのかもしれないけど、それはそれで確かめる必要があると思っている。

めくるめく符号と暗号の世界・・・ｗ


# 20230113

昨日もよく眠れなかった。

3時間しか寝てないので集中できず、思わずGF(p^m)などの一般化した拡大帯の演算について実装してました。

ここで予期していなかったのは、2の拡大体の場合はXORだけで完結していたのが、一般の足し算の場合は、その拡大体を構成する原始多項式を法とするような演算になるということに気がついたことです。

そしてそれは半年くらい前に実装していたvan-koran.cに凝縮されています。

一般の拡大体が余り応用で見られないのは、演算が高速化できないからではないかと思い、ここでなにかいいアイデアがあればそれも成果になりそうです。
もし2次元符号の研究が失敗したら、スターもついたことだし、あえて誰も取り組みそうにないワイルドな方のこの研究を進めてみたいと思います。

とはいえボケた頭では何も理解できないのですがｗ

# 20230111

ほとんど寝ないで作業。
と行っても真夜中に目が覚めただけなんだけど。
確か昔やったような計算が合ったなと思いだして実際動かしてみる。

５の５次拡大退場の計算では、整数をそのままかけたり足したりしてもうまく行かない。
これらの演算は、GF(5)上の多項式、つまり係数が５までの多項式の四則演算として扱わなければならない。

ここでまた私が作った手前味噌の多項式らぶらりが役に立つ。
このライブラリは自分で作ったものとは思えないほど大活躍している。
超楕円曲線の離散対数も、この手前味噌のライブラリが大活躍したのだ。

それを思い出して早速計算したら正解だった。
そしてBM法で５の拡大体の誤り訂正を試して見る予定。（パターソンは面倒なのでやらない）

最初はもうダメかと諦めていたけど、多項式でなく整数の演算ばかりやっていたのだからうまく行くはずがない。

そんな感じで、言うわけであとは寝て起きたら続きをします。

# 20230109

いやー、もう何年やってるんだかという感じですが、クラシカルGoppaは行き詰まって疲れたので別のテーマをかじってみることに。
Wild Niederreiterにはパターソン復号の素晴らしい一般化があるので、今までの延長でできそうな気がするし、それがダメなら群符号か
電子署名がやってみたい。

特に群符号が未知の領域であり、ほんと今まで聞いたことがあってなんとなく気になるから論文だけ持っているという状態。

息抜きに群符号でちょっと遊んでみるのもいいかもしれない。
２次元符号はグレブナ入ってるし、グレブナなしの群符号ってないの？
という意味でもちょっと読もうとしている。

署名は攻撃法と合わせて考えないといけないのでちょっと重そうな感じがして気が引ける。
気のせいかもしれないけど。

なのでちゃんと目を通してないからわからないけど、ワイルドか、署名か、群符号のどれかにする予定！

### 追記：
ちょっと読んでみました。
群符号とはいっても別に普通のBCHみたいなのをいうだけであって、他の軍符号は全部リードソロモンとBCHに取って代わられてしまい、
普通の古い富豪になってしまっている。
私の誤った認識では、群符号というのは２次元符号のことだったのだ。
しかし、これを復号しようとすると、既存の論文ではすべてグレブナを必要とする。
ちょっと前まで２変数の割り算をやっていたので、その勢いでグレブナに突入できればいいけど、今はちょっとそこまでの勢いがない。
グレブナができるようになれば、２次元符号も代数幾何符号も怖くないんですが。

で、ワイルドも読んでみたけど結構難しかった。

弱ポポフ形式って何？ｗ

既存の知識に
$+ \alpha$
って感じでしょうかｗ

これがちょっと重荷で辛い。

そして最後の切り札として用意されていたのが、暗号学的プロトコル。
ゼロ知識証明とかです。
ゼロ知識証明はもう暗号やってる人には半ば常識のようなところがあるので、もっと珍しい紛失通信や、さらに必殺技をやるかもしれない。
でもなんだか興味が持てない。
大人しく４でもやろうかな。

まあ、迷いがあるときは調子の悪いときでもある。

ダブルサイズのシンドロームへの情熱を、グレブナや弱ポポフ形式にかけることができるかどうか、全てはそこにかかっている。
弱ポポフ形式を理解するより、多項式の割り算のほうが簡単なような気がするのだが。

証拠識別不可能性とか、これも読んでみないとわからないですが、とりあえず明日以降決めます。
あともう一つはひみつですｗ

# 20221028

いやー、もう年末ですね。

年明けから狙っていたはずのWild-Goppaですが、精神状態が安定しないので来年に持ち越すことになりそうです。
というより、暗号の何がそんなに面白かったのかわからなくなってきたという感じです。
これも病気のせいだと思うんですが、また思い出したら突発的に再開するかもしれないです。
いまルバートさんというC言語の先生に色々お手本を見てもらいながら、リファクタリングする予定です。
このままコードを公開しているだけでは、誰にも使ってもらえないし、それではなんのためのコードなのか意味がないからです。
今は代数幾何的符号に興味がありますが、ゴールの近さから言うと符号ベースの電子署名といったところです。

このテーマを続けていく上で困難なのは、基礎体上の計算でBM法やパターソン法が、ほぼそのまま使えるのかどうかです。
そこを乗り越えていかないとこの課題はこなせない。
しかしまだその見通しはついていかない感じです。

その代わり代数幾何符号の部分対部分符号はMcEliece暗号で使える新たな符号クラスであることが比較的最近わかっており、
今はその基礎である一般代数曲線の1点生成符号の復号法に興味を持っているところです。
まあ、学部生の卒検レベルに相当するかしないかという感じです。

お金がないので来月まで引きこもっていようと思うのですが、数式が頭に入ってこない状態なのでソーカル事件の本とかスノーデンの本とかに目を通しています。
自分たちは、資本主義という経済システムの中で決まる価値観を、内面化しているのだという。
その上で、人間の命の価値までが労働生産性のみで決まってしまうということに、生命に対する侮辱を感じてしまうのであった。

分析哲学面白そうｗ

# 20220902

どこを書き足したかよくわからなかったので追記。  
1331の元を、曲線ha[5][4]で有理点をカウント。  
壊れてません。
計算が終わるのにAMD2700Xで２分弱かかります。  
カウントの数は30589が正解。  
odd.cがメインディッシュ。

# 20220821

実際に動かしてみてわかったのは、足し算が複雑になって遅いという事でした。

足し算ならバイナリだと論理演算が１つで何もすることがないのですが、３以上の基礎体を扱う場合はそれらの原始多項式を法とする多項式演算が必要で、
あらかじめ演算テーブルとして持つことでなんとかしてます。
実行するごとに作ると無駄に待ち時間がかかるので、多項式と整数の対応表を作ろうと思います。

とりあえずGF(3^3)で代数曲線の点を数えることができるようになりました。

# 20220819

普通にGoppaやってればよかったのに、AGにまで無駄に時間を費やしてしまった。
11の拡大帯に変更するのって、足し算関数が必要なので結構面倒だし計算速度も遅い。

# 20220816

バイナリAG符号、または奇標数の拡大体でのAG符号（代数幾何符号）は攻撃法が見つかっていないので未だに安全だという論文を見つけた。
符号ベースの暗号に使える符号が増えたのは誠に喜ばしい。
古典的な方も素体上での攻勢は知らないのでちょっとやってみるつもり。

# 20220815

暗号・・・を作りた買ったのですが、全然形になってないですし終わる見通しも立ちません。


# 20220731

凄い熱気。外はまるでサウナのよう。調子が悪くて全然先に進めませんでした。

