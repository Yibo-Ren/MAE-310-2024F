 R = 0.3; //圆的半径
L = 1.0; //矩形边长一半

Point(1) = {L, -L, 0}; 				       //点（1）位置，对应在板的右下角
Point(2) = {L, L, 0};				       //点（2）位置，对应在板的右上角
Point(3) = {-L, L, 0}; 				       //点（3）位置，对应在板的左上角
Point(4) = {-L, -L, 0};                                //点（4）位置，对应在板的左下角的圆心
Point(5) = {-L + R, -L, 0};			       //点（5）位置，对应在板底边与圆弧交点
Point(6) = {-L, -L + R, 0};			       //点（6）位置，对应在板左边与圆弧的交点
Point(7) = {-L + Cos(Pi/4) * R, -L + Sin(Pi/4) * R, 0};//点（7）位置，对应在板对角线与圆弧的交点

Circle(1) = {5, 4, 7};//在点5，4，7之间生成圆弧，编号1
Circle(2) = {7, 4, 6};//在点7，4，6之间生成圆弧，编号2

Line(3) = {6, 3};//在点6，3之间生成一条直线，编号3
Line(4) = {3, 2};//在点3，2之间生成一条直线，编号4
Line(5) = {2, 1};//在点2，1之间生成一条直线，编号5
Line(6) = {1, 5};//在点1，5之间生成一条直线，编号6
Line(7) = {2, 7};//在点2，7之间生成一条直线，编号7


//+
Curve Loop(1) = {4, 7, 2, 3};//由曲线4，7，2，3构成封闭曲线
Plane Surface(1) = {1};	     //由封闭曲线构成平面

Curve Loop(2) = {7, -1, -6, -5};//由曲线7，-1，-6，-5构成封闭曲线，有方向，负号代表由大编号点指向小编号点
Plane Surface(2) = {2};		//由封闭曲线构成平面


Transfinite Line{1, 2, 3, 4, 5, 6, 7} = 3; //对边划分网格

Transfinite Surface{1};			   //对面划分网格
Transfinite Surface{2};			   //对面划分网格

Recombine Surface{1};			   //重组为四边形网格
Recombine Surface{2};

Mesh.ElementOrder = 1;			//网格元素阶数为1阶线性
Mesh.Algorithm = 8;			//网格细化

// EOF
//+


//+
Physical Curve("x_axis", 8) = {6};
//+
Physical Curve("y_axis", 9) = {3};
//+
Physical Curve("wall", 10) = {4};
//+
Physical Curve("boundary", 11) = {5};
//+
Physical Curve("hole", 12) = {2, 1};
