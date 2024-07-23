from typing import Callable, Optional, Any
import asyncio

import cv2, zxingcpp

import numpy as np

from nicegui import ui

import base64


class Camera():
    def new__init__(
            self,
            icon: str,
            icon_color: str,
            background_color: str,
            for_id: Optional[str] = 'fileupload',
            canvas_id: Optional[str] = 'canvas',
            on_change: Optional[Callable[..., Any]] = None,
            compression: Optional[float] = 0.8,
    ):
        self.compression = compression
        self.canvas_id = canvas_id
        self.for_id = for_id
        self.icon = icon
        self.icon_color = icon_color
        self.background_color = background_color
        self.on_change = on_change

        ui.add_body_html(f'''
        <canvas id="{self.canvas_id}" style="display:none;"></canvas>
        <script>
        function setupCamera{self.for_id}() {{
            const fileInput = document.getElementById('{self.for_id}');
            const label = document.querySelector('label[for="{self.for_id}"]');

            if (fileInput && label) {{
                label.addEventListener('click', function(event) {{
                    event.preventDefault();
                    fileInput.click();
                }});

                fileInput.addEventListener('change', handleFileSelect_{self.for_id});
            }}
        }}

        function handleFileSelect_{self.for_id}(event) {{
            const file = event.target.files[0];
            if (!file || !file.type.match('image.*')) {{
                alert('Please select an image file.');
                return;
            }}

            const reader = new FileReader();
            reader.onload = function(readerEvent) {{
                const img = new Image();
                img.onload = function() {{
                    const maxWidth = 2048;
                    const maxHeight = 1536;

                    let width = img.width;
                    let height = img.height;
                    if (width > height) {{
                        if (width > maxWidth) {{
                            height *= maxWidth / width;
                            width = maxWidth;
                        }}
                    }} else {{
                        if (height > maxHeight) {{
                            width *= maxHeight / height;
                            height = maxHeight;
                        }}
                    }}
                    const canvas = document.getElementById('{self.canvas_id}');
                    const ctx = canvas.getContext('2d');
                    canvas.width = width;
                    canvas.height = height;
                    ctx.clearRect(0, 0, canvas.width, canvas.height);  // Clear the canvas
                    ctx.drawImage(img, 0, 0, width, height);

                    // Call the function to process the canvas image
                    processImage_{self.for_id}();

                    // Reset the file input
                    event.target.value = "";
                }};
                img.src = readerEvent.target.result;
            }};
            reader.readAsDataURL(file);
        }}

        function processImage_{self.for_id}() {{
            const canvas = document.getElementById('{self.canvas_id}');
            const imageURL = canvas.toDataURL("image/jpeg", {self.compression});
            // Call Python callback
            pywebview.api.on_image_processed('{self.for_id}', imageURL);
        }}

        // Setup the camera when the DOM is fully loaded
        if (document.readyState === 'loading') {{
            document.addEventListener('DOMContentLoaded', setupCamera{self.for_id});
        }} else {{
            setupCamera{self.for_id}();
        }}
        </script>
        ''')

        self._camera_input = ui.input(on_change=self.on_change)
        self._camera_input.props(f'type="file" accept="image/*" capture="camera" id="{self.for_id}"')
        self._camera_input.classes('w-full').style('display: none;')

    def __init__(
            self,
            icon: str,
            icon_color: str,
            background_color: str,
            for_id: Optional[str] = 'file-upload',
            canvas_id: Optional[str] = 'canvas',
            on_change: Optional[Callable[..., Any]] = None,
            compression: Optional[float] = 0.5,
    ):
        """Camera Input
                :param icon: The icon to be used for the camera
                :param icon_color: The color of the icon
                :param background_color: The color of the background, must be in rgba format
                :param for_id: The id of the input being used, this should be set if you're using this more than once on a page (default: file-upload)
                :param canvas_id: The id for the canvas, this should be set if you're using more than 1 canvas on the page (default: canvas)
                :param on_change: callback to execute when the value changes
                :param compression: Compression ammount used on uploaded images, change this if no image is being returned usually caused by large file sizes (default: 0.8)
        """
        self.compression = compression
        self.canvas_id = canvas_id
        self.for_id = for_id
        self.icon = icon
        self.icon_color = icon_color
        self.background_color = background_color
        ui.add_body_html(f'''
        <canvas id="{self.canvas_id}" style="display:none;"></canvas>
        <script>
        window.addEventListener("DOMContentLoaded", function() {{
            document.getElementById('{self.for_id}').addEventListener('change', function(event) {{
                const file = event.target.files[0];
                if (!file.type.match('image.*')) {{
                    alert('Please select an image file.');
                    return;
                }}
                
                console.log(file);

                const reader = new FileReader();
                reader.onload = function(readerEvent) {{
                    const img = new Image();
                    console.log(img);
                    img.onload = function() {{
                        const maxWidth = 2048;
                        const maxHeight = 1536;

                        let width = img.width;
                        let height = img.height;
                        if (width > height) {{
                            if (width > maxWidth) {{
                                height *= maxWidth / width;
                                width = maxWidth;
                            }}
                        }} else {{
                            if (height > maxHeight) {{
                                width *= maxHeight / height;
                                height = maxHeight;
                            }}
                        }}
                        const canvas = document.getElementById('{canvas_id}');
                        canvas.width = width;
                        canvas.height = height;
                        const ctx = canvas.getContext('2d');
                        ctx.drawImage(img, 0, 0, width, height);
                    }};
                    img.src = readerEvent.target.result;
                }};
                reader.readAsDataURL(file);
            }});
        }}, false);
        </script>
        ''')
        self._camera_input = ui.input().on("update:model-value", on_change) #(on_input=on_change)
        self._camera_input.props(f'type="file" accept="image/*" capture="camera" for="{self.for_id}"')
        self._camera_input.classes('w-full').style('display: none;')

    def show_camera(self):
        ui.html(rf'''
        <label for="{self.for_id}" class="custom-{self.for_id} q-icon text-{self.icon_color} notranslate material-icons q-pa-sm" style="cursor: pointer; font-size: 64px; background-color: {self.background_color}; border-radius: 10px;" aria-hidden="true" role="presentation" rounded="true">
            {self.icon}
        </label>
        ''')

    def _on_change(self):
        if self.on_change:
            self.on_change()

    async def get_image(
            self,
            on_change: Optional[Callable[..., Any]] = None,
    ):
        await asyncio.sleep(0.5)
        file = await ui.run_javascript(
            f'var canvas = document.getElementById("{self.canvas_id}"); var imageURL = canvas.toDataURL("image/jpeg", {self.compression}); return imageURL;',
            timeout=15)
        return file
