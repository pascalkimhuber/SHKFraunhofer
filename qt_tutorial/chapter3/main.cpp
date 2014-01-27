#include <QApplication>
#include <QPushButton>
#include <QFont>
#include <QWidget>
#include <QVBoxLayout>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	
	QWidget box;
	box.resize(10, 12);

	QPushButton quit("Quit", &box);
	quit.setFont(QFont("Times", 100, QFont::Bold));

	QObject::connect(&quit, SIGNAL(clicked()), &a, SLOT(quit()));
	
	box.show();

	return a.exec();
}










