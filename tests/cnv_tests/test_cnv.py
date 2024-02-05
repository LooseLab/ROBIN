from methnicegui.copy_number_component import index_page
from nicegui.testing import Screen


def test_label_message(screen: Screen) -> None:
    """
    This tests that the cnv element is present on the page.
    """
    index_page()
    screen.open('/')
    screen.should_contain('Copy Number Variation')

def test_CNV_Plot_label(screen: Screen) -> None:
    """
    This tests that the BAM file is loaded and the correct bin width is calculated and displayed.
    """
    index_page()
    screen.open('/')
    screen.should_contain('Bin Width: 19546116')

def test_button_click(screen: Screen) -> None:
    """
    Generic test for the functionality of the button in the footer component.
    """
    index_page()
    screen.open('/')
    screen.click('More Information')
    screen.should_contain('Looselab')
